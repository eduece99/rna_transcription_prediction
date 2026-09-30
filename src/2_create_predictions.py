


import pandas as pd
import numpy as np
import re
import xgboost as xgb
from pathlib import Path

from utils import generate_features_from_dataframe


def make_predictions(model_med : xgb.Booster, model_low : xgb.Booster, model_high : xgb.Booster, X : pd.DataFrame) -> dict[str, np.ndarray]:
    """Generate model predictions for a feature table.

    Args:
        model_med: Trained XGBoost regressor for the median outcome.
        model_low: Trained low-quantile XGBoost model.
        model_high: Trained high-quantile XGBoost model.
        X: DataFrame of input features for prediction.

    Returns:
        Dictionary with index, lower80, upper80, and mean_liver_rna predictions.
    """
    
    predictions = {
        "index": X.index,
        "lower80": model_low.predict(X),
        "upper80": model_high.predict(X),
        "mean_liver_rna": model_med.predict(X),
    }
    return predictions




if __name__ == "__main__":

    utr_df = pd.read_csv("utr_library.csv")
    invivo_df = pd.read_csv("invivo_liver.csv")

    feature_df = generate_features_from_dataframe(utr_df)
    df_features = pd.concat([utr_df, feature_df], axis=1)

    utr_liver = pd.merge(utr_df, invivo_df, on="construct_id", how="inner")

    # no model should be trained on means, but may be useful as a best case accuracy guide
    utr_liver_means = utr_liver.groupby(["construct_id", "day"]).agg(
        mean_liver_rna=("liver_rna", "mean"),
        lower20=("liver_rna", lambda x: np.percentile(x, 20)),
        upper80=("liver_rna", lambda x: np.percentile(x, 80))
    ).reset_index()
    # somehow, I think that using day 3 in addition to day 28 might be useful, but let's start with just day 28 for now

    # features of relevant constructs
    
    day28_features_means = pd.merge(
        utr_liver_means[utr_liver_means["day"] == 28],
        df_features,
        on="construct_id",
        how="inner"
    )



    # min max scaling of data?

    # cohort and day are the same and thus redundant.  Cage Id however, is not

    X_means, y_means = day28_features_means[feature_df.columns].set_index(day28_features_means["construct_id"]  ), day28_features_means["mean_liver_rna"]
   
    model_med = xgb.XGBRegressor()
    model_med.load_model("models/model_med.json")
    model_low = xgb.XGBRegressor()
    model_low.load_model("models/model_low.json")
    model_high = xgb.XGBRegressor()
    model_high.load_model("models/model_high.json")

    predictions = pd.DataFrame( make_predictions(model_med, model_low, model_high, X_means) ).set_index("index")
    predictions.to_csv("predictions.csv")


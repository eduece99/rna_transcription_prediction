"""

want liver rna levels after 28 days - invivo_liver.csv prediction
sequence info from prediction_targets.csv/utr_library.csv

Typically 5 mice per screen, so average in vivo liver RNA levels across 5 mice for each construct.
(though best to do a check on this and filter out exceptions - missing data)

Also report 20th and 80th percentile of the distribution of liver RNA levels across the 5 mice for each construct.


"""


import pandas as pd
import numpy as np
import re
import xgboost as xgb
from sklearn.model_selection import GroupKFold, cross_val_score
from sklearn.model_selection import RepeatedKFold
from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error, median_absolute_error
from pathlib import Path
#import typer

from utils import generate_features_from_dataframe

# PROMPT  please provide me some code which splits the RNA sequences into features





# Inspired from https://scikit-learn.org/stable/auto_examples/ensemble/plot_gradient_boosting_quantile.html
# also https://xgboost.readthedocs.io/en/stable/python/examples/prediction_intervals.html#sphx-glr-python-examples-prediction-intervals-py
# AI responses seemed to convoluted so I did this manually with respect to the 2 references
def fit_models_liver_rna(X : pd.DataFrame, y: pd.Series, groups: pd.Series) -> pd.DataFrame:
    """Train median and quantile XGBoost models for day-28 liver RNA.

    Args:
        X: Feature matrix for construct-level or mouse-level predictors.
        y: Target vector of liver RNA readouts.
        groups: Group labels for group-based cross-validation.  Needed if multiple measurements per construct (e.g., multiple mice per construct).

    Returns:
        Tuple containing the median model, lower-quantile model, and
        upper-quantile model.
    """

    common_params = {
        "n_estimators": 100,
        "max_depth": 5,
        "learning_rate": 0.2,
        "subsample": 0.8,
        "colsample_bytree": 0.8,
        "random_state": 42
    }

    model_med = xgb.XGBRegressor(
        objective="reg:squarederror",
        **common_params
    )

    model_low = xgb.XGBRegressor(
        objective="reg:quantileerror",
        quantile_alpha=0.10,
        **common_params
    )

    model_high = xgb.XGBRegressor(
        objective="reg:quantileerror",
        quantile_alpha=0.90,
        **common_params
    )

    model_med.fit(X, y)
    model_low.fit(X, y)
    model_high.fit(X, y)

    # report statistics - performance etc

    

    
    # OOF predictions
    gkf = GroupKFold(n_splits=5)
    groups = day28_features["construct_id"]

    pred_med = np.zeros(len(X))
    pred_low = np.zeros(len(X))
    pred_high = np.zeros(len(X))
    fold_r2 = []

    for train_idx, val_idx in gkf.split(X, y, groups):
        X_tr, X_val = X.iloc[train_idx], X.iloc[val_idx]
        y_tr, y_val = y.iloc[train_idx], y.iloc[val_idx]

        model_med.fit(X_tr, y_tr)
        model_low.fit(X_tr, y_tr)
        model_high.fit(X_tr, y_tr)

        pred_med[val_idx] = model_med.predict(X_val)
        pred_low[val_idx] = model_low.predict(X_val)
        pred_high[val_idx] = model_high.predict(X_val)
        fold_r2.append(r2_score(y.iloc[val_idx], pred_med[val_idx]))

    print("Outputting prediction statistics for day 28 liver RNA levels from XGBoost models:")

    print(f"Fold R2 values: {[round(score, 3) for score in fold_r2]}")

    prediction_statistics(y, pred_med, pred_low, pred_high)

    return model_med, model_low, model_high


def prediction_statistics(y_true: pd.Series, pred_med: pd.Series, pred_low: pd.Series, pred_high: pd.Series) -> dict[str, float]:
    """Compute baseline prediction statistics.

    Args:
        y_true: True target values.
        pred_med: Median predicted target values.
        pred_low: Lower quantile predicted target values.
        pred_high: Upper quantile predicted target values.

    Returns:
        Dictionary containing R2, MAE, RMSE, and mean residual.
    """

    
    residual = y_true.to_numpy() - pred_med
    abs_error = np.abs(residual)
    q1, q3 = np.quantile(residual, [0.25, 0.75])
    iqr = q3 - q1
    outlier_mask = (residual < q1 - 1.5 * iqr) | (residual > q3 + 1.5 * iqr)


    print(f"OOF R2: {r2_score(y_true, pred_med):.3f}")
    print(f"MAE: {mean_absolute_error(y_true, pred_med):.3f}")
    print(f"RMSE: {np.sqrt(mean_squared_error(y_true, pred_med)):.3f}")
    print(f"Mean residual (observed - predicted): {residual.mean():.3f}")
    print(f"Mean 80% interval width: {np.mean(pred_high - pred_low):.3f}")
    

    diagnostics = pd.DataFrame(
        {
            "observed": y_true.to_numpy(),
            "predicted": pred_med,
            "residual": residual,
            "absolute_error": abs_error,
            "iqr_outlier": outlier_mask,
        },
        index=y_true.index,
    )
    print("Largest errors / IQR outliers:")
    print(diagnostics.loc[outlier_mask].sort_values("absolute_error", ascending=False))

    #return stats


if __name__ == "__main__":
    
    utr_df = pd.read_csv("utr_library.csv")
    invivo_df = pd.read_csv("invivo_liver.csv")

    feature_df = generate_features_from_dataframe(utr_df)
    df_features = pd.concat([utr_df, feature_df], axis=1)

    utr_liver = pd.merge(utr_df, invivo_df, on="construct_id", how="inner")

    # features of relevant constructs
    day28_features = pd.merge(
        utr_liver[utr_liver["day"] == 28],
        df_features,
        on="construct_id",
        how="inner"
    )


    # means for baseline purposes - not used for model training
    group_keys = ["construct_id", "day"]
    grouped_rna = utr_liver.groupby(group_keys)["liver_rna"]

    utr_liver_stats = utr_liver.assign(
        mean_liver_rna=grouped_rna.transform("mean"),
        lower20=grouped_rna.transform(lambda values: values.quantile(0.20)),
        upper80=grouped_rna.transform(lambda values: values.quantile(0.80)),
    )
    # somehow, I think that using day 3 in addition to day 28 might be useful, but let's start with just day 28 for now

    # features of relevant constructs
    
    day28_features_means = pd.merge(
        utr_liver_stats[utr_liver_stats["day"] == 28],
        df_features,
        on="construct_id",
        how="inner"
    )

    

    #X_means, y_means = day28_features_means[feature_df.columns].set_index(day28_features_means["construct_id"]  ), day28_features_means["mean_liver_rna"]
    X, y = day28_features[feature_df.columns], day28_features["liver_rna"]

    print("Outputting prediction statistics for day 28 liver RNA levels against the means of day28 data:")
    prediction_statistics(y, day28_features_means["mean_liver_rna"], day28_features_means["lower20"], day28_features_means["upper80"])
    
    

    #model_med.score(X, y)  # R^2 score on training data

    if not Path("models").exists():
        Path("models").mkdir(parents=True)

    model_med, model_low, model_high = fit_models_liver_rna(X, y, day28_features["construct_id"] )

    model_med.save_model("models/model_med.json")
    model_low.save_model("models/model_low.json")
    model_high.save_model("models/model_high.json")
 


    

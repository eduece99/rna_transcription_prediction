# Report

## Usage

To deploy the python environment (installation), please run `./deploy.sh` (or `bash deploy.sh` if not executable).  

Then, `./run_scripts.sh` to run the Python scripts to output the predictions into predictions.csv, and the suggested new constructs into designs.csv

## Recommendation and interpretation

10 custom constructs are in the "designs.csv" file with predicted RNA levels and 80% lower and upper prediction intervals.  The same predicted levels have been outputted for the existing constructs to act
as a baseline.  The baseline here is of course a fair comparison, as the data are already there and have been experimentally validated.  

Said constructs maintain a strong Kozak-like 5'UTR context while avoiding obvious decay liabilities in the 3'UTR. In this setting, the key decision is not only which sequences predict the highest day-28 RNA, but which constructs are most likely to stay durable after a month in liver without being crippled by AU-rich decay elements or the liver-dominant miR-122 seed motif.

The workflow in the three analysis scripts is intentionally simple and transparent: first build sequence-derived features from the 5'UTR and 3'UTR, then fit a mean XGBoost regressor and quantile XGBoost models for uncertainty, and finally use a mutation-based genetic search to propose new candidates under the learned model. This matches the biology described in the brief: 5'UTR effects are dominated by translation and initiation, while 3'UTR effects are dominated by transcript stability and decay.

## Methodology

The modeling pipeline is split across three scripts:

- 1_create_models.py trains the predictive models. It reads the UTR library and mouse liver data, builds sequence features for each construct, and fits an XGBoost model for day-28 liver RNA. The mean model is used for the expected value, and quantile models generate lower and upper bounds for an uncertainty band.
- 2_create_predictions.py loads the trained models and predicts day-28 mean values for the target constructs using the same engineered feature set. This provides the predicted mean and approximate 80% interval for each construct.
- 3_create_constructs.py uses a mutation-based search over valid UTR sequences, scoring each new candidate with the trained median model and heavily penalizing biologically poor motifs such as upstream AUGs in the 5'UTR, AU-rich decay motifs like AUUUA, and the miR-122 seed ACACUCC in the 3'UTR.



The feature engineering focuses on the features explicitly highlighted in the brief: upstream AUG counts and Kozak context in the 5'UTR; AU-rich motifs, miR-122 seed matches, and composition-based decay proxies in the 3'UTR. This is not an opaque end-to-end model — it is a biologically informed predictor built around the mechanisms the exercise is designed to reward.

## Rationale

XGBoost has been a well known library for the highly performant gradient boosted tree method, which can be used to yield an excellent performance on a dataset without requiring a great training time.  Granted it is not as transparent as say, a linear regression model, but interpretability wasn't my main concern here.  In addition, feature importances, if needed, can be extracted directly from the tree model anyway by virtue of how soon and thus influential (smaller depth number) a feature is used in a tree, collectively.

More expensive methods such as neural networks, I avoided, namely due to hardware constraints.  I am aware of more recent tabular predictive methods such as TabPFN, but given their compute time I did not see any real reason to use them for a small predictive gain.  In addition, the data used here was somewhat lacking and it would be more useful to obtain better data, as explained later.

I considered a plain “best fit” model alone, but that would give a point estimate without an honest uncertainty band. I also considered a purely random sequence generator, which would not be constrained by the biology we know matters here. The final workflow is a fair compromise: the model is read out in the same biological language as the domain problem, and the design-generation step is constrained by sequence validity and known negative motifs.

A basic baseline would be a simple mean model or a library-only heuristic such as “avoid AUUUA and ACACUCC.” That is transparent and easy to explain, but it ignores the richer combination of 5'UTR initiation effects and 3'UTR decay effects. The final XGBoost approach is more informative because it learns how these features jointly contribute, while still remaining interpretable enough to explain in a bench-scientist setting.

## How much confidence I have

Statistics on latest computation of models:

OOF R2: 0.394
MAE: 1768.501
RMSE: 3283.903
Mean residual (observed - predicted): -249.415
Mean 80% interval width: 2444.588
Fold R2 values: [0.535, -1.615, -1.448, 0.714, 0.468]

Compared however, to just taking the means and using these as predictions for each construct of the in vivo data, as a baseline:
OOF R2: 0.841
MAE: 758.726
RMSE: 1681.263
Mean residual (observed - predicted): 0.000

To summarise, they're predictive, but weak compared to just taking means.  I did not have enough time to figure out if the differences between the two "methods" are significant.  Of course, you cannot use the means alone as predictions, as they are based on already tested data.

The intervals are based on quantile XGBoost models trained on the same feature set and then interpreted in the context of the actual decision: the mean of a new five-mouse study using freshly prepared material. 

That said, the uncertainty still reflects a limited dataset. I would treat the interval as a directional decision aid, especially for constructs near the top of the ranking. The model should be trusted most for relative ranking and for identifying poor motif patterns, not for claiming exact absolute values on a novel construct.

## What would change my recommendation

In short, more data and more time.

5 mice per construct for in vivo data is useful, though considering the amount of variability, more would always help to establish a better mean prediction.  I could recommend 30 as a rule of thumb for normally distributed data, though of course that is expensive.

I did not use the in vitro data due to time constraints.  However more crucially, it was unclear whether these RNA levels were measured after 28 days, which was a specified prerequisite of the predictions.  Thus, I ignored these data.

The modelling was performed only using day 28 data points, ignoring day 3 from the in vivo data.  If say, the ratio between day 3 vs day 28 was taken into account (and thus model degradation) then this would provide more useful predictions for future mRNA.

Finally, The prediction intervals are around the median, not the mean, owing to being models trained on estimating quantiles.  This was due to a time constraint - I didn't have time to properly introduce
a bootstrapping based estimation for estimating intervals around the mean, which would have been a more informative prediction interval.  This means that on occasions, the mean prediction can actually be outside one of the prediction intervals.


## What I would measure next

In vitro data after days 3 and 28 on the same constructs as used in the in vivo data at a minimum, as an independent experiment.  This would reduce potential noise from variation in mice, plus the day 3 vs day 28 ratio would give me a useful measure of degradation, which I have not directly modelled.  

If the in vitro and in vivo data generally agree, then more constructs can be tested in vitro.  As mentioned earlier, given the ambiguity on which day the in vitro data was after, I did not use it.
 



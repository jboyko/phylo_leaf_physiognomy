# Climate models

All prediction traits come from the site being predicted. Leaf measurements are averaged within species × site; fossil traits and ages are never pooled across sites. Training can still use species means across extant calibration sites.

| Model | Fitted to | Site estimate | Current fossil application |
| --- | --- | --- | --- |
| LM: species | Species means | Predict local species → average | Validation only |
| PGLS | Species means, accounting for phylogeny | Predict local species → average | Validation only |
| PIP | Same fit as PGLS | Predict local species + phylogenetic correction → average | Primary reconstruction |
| 12-trait site regression | Site means including appropriate untoothed tooth values | Average local traits → predict | Fossil baseline |
| LM: site, Peppe-style | Site means excluding missing untoothed tooth measurements | Average local traits → predict | Validation only |
| DiLP regression | Three published predictors per target; coefficients refitted within folds | Published site averaging → fitted equation | Validation only |

The five original fitted approaches have imputed-training and complete-case-training validation variants. PGLS/PIP impute missing prediction inputs in both variants. Current fossil predictions use imputation learned from extant calibration data. Biological tooth-value filling is separate from statistical imputation.

MAT predictions are averaged arithmetically; MAP predictions are averaged on the log scale and exponentiated, giving a geometric mean in cm.

`specimen` and `sp_zero` are duplicate site-aggregation configurations in the current code. Fossil outputs retain `lm_site` and `pip_site`. Cross-site pooled fossil comparators (`lm_sp`, `pip_sp`) and their placement records were removed; species-trained validation models remain because their held-out predictions already use local traits.

DiLP regression uses complete finite inputs. The Dana comparison scores all three headline models on the same eligible sites within each climate target. Run `code/03c_dilp_cv.R` after the main CV to reproduce these matching-site scores. The climate calibration contains 92 sites, with Yasuni ridgetop and upper slope combined before averaging traits and assigning folds.

Cross-validation RMSE measures error across held-out extant sites. Site-specific fossil prediction intervals require the joint covariance of species prediction errors; see [prediction uncertainty](prediction_uncertainty.md).

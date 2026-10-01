# MAP model experiment — 24 September 2026

## Result

Nested recalibration of PIP was the best candidate in this exploratory comparison:
log(MAP) RMSE fell from 0.5295 to 0.4220 (20.3% lower RMSE; 36.5% lower MSE).
It improved RMSE in eight of ten outer folds. Site-trained random forest gave a
smaller improvement (0.4963); the tested GAMs and species-trained RF did worse
than PIP. These results concern point prediction, not prediction-interval coverage.

| Model | RMSE, log cm | MAE | Mean error | Predicted / observed SD |
|---|---:|---:|---:|---:|
| PIP + recalibration | 0.4220 | 0.3151 | 0.0161 | 0.707 |
| Site RF | 0.4963 | 0.3565 | 0.0095 | 0.591 |
| PIP + offset only | 0.5051 | 0.3898 | 0.0070 | 0.318 |
| PIP | 0.5295 | 0.4129 | 0.1576 | 0.318 |
| Species RF | 0.5867 | 0.4628 | 0.1832 | 0.286 |
| Species GAM | 0.6022 | 0.4815 | 0.2007 | 0.261 |
| Site GAM | 0.6060 | 0.4308 | 0.0138 | 0.665 |
| Species LM | 0.6270 | 0.4989 | 0.2170 | 0.228 |
| Mean only | 0.6464 | 0.5216 | 0.0002 | 0.034 |
| Site LM | 0.6879 | 0.4584 | 0.0430 | 0.815 |

## What recalibration changes

After averaging the species log(MAP) predictions within each site, the correction
is `corrected_log_MAP = a + b * predicted_log_MAP`. Within each outer training
set, five inner site-grouped folds supply out-of-fold PIP site predictions.
Regressing their observed log(MAP) on those predictions learns a and b; these
coefficients then transform the outer held-out PIP site predictions.

The ten slopes range from 1.95 to 2.48 (mean 2.18). Prediction SD increases from
0.204 to 0.455 log units, versus observed SD 0.643. Correlation changes little
(0.750 to 0.753), consistent with correcting the scale of an existing signal.
An offset-only control, estimated from the same inner predictions, reaches RMSE
0.5051. Thus correcting the mean bias alone does not explain the full gain.
These observations do not identify the biological or statistical cause of the
original compression.

## Design

- Existing ten site-grouped outer folds, 92 sites total, with roughly nine sites
  excluded at a time. This benchmark PIP RMSE is 0.5295; the separately run
  leave-one-site-out benchmark was 0.5327. Every method here uses the same ten folds.
- Species-level and site-level trait-only LM, GAM, and random forest use the
  existing fossil-measurable predictors, tooth filling, and aggregation rules.
  Species predictions are averaged on the log scale; site-trained models use
  local site means. No cross-site traits enter a held-out site's prediction.
- Bagged imputation, active-predictor filtering, scaling, and constant-predictor
  filtering are fitted on each respective training subset. Species-level
  calibration retains the existing operational taxon labels and scaffold filter.
- PIP retains the production phylogenetic correction and scaffold covariance.
  New GAM/RF candidates are trait-only; no nonlinear phylogenetic model was tested.
- GAM: additive shrinkage smooths (`bs="ts"`), k=3 or k=4 per nonconstant trait,
  REML smoothing with gamma=1.4; inner site MSE chooses k.
- RF: `ranger` 0.18.0; mtry=3 or 8, minimum node size=3 or 10; 400 trees for inner
  evaluation and 800 for outer fits; inner site MSE chooses the combination.
- Five randomly balanced inner site folds for each outer training set, with
  deterministic recorded seeds. Fifty inner PIP fits learn recalibration; outer
  PIP uses the saved matching production fits. Inner tuning losses weight sites
  equally. The offset-only control also uses only inner held-out residuals.

## Verification

All 92 outer predictions are finite. Audits verify disjoint outer and inner site
memberships, no outer test site anywhere in its inner folds, and exactly one
inner held-out prediction per training site per candidate. Recomputing selected
settings and recalibration from saved inner predictions reproduces the output.
PIP and both LM controls reproduce the existing benchmark predictions within
1e-7 (LM differences were approximately 1e-14). Observed climate is checked
against `data/dat_site.csv`. The result script refuses incomplete checkpoints.

The run emitted the existing input-data warnings about outlier measurements and
nonbinary margin states; no model convergence or prediction failures occurred.
`dilp` applies its existing conversion of nonbinary margin values to 0.5.

## Interpretation and limits

Recalibration merits follow-up before spending more effort on model flexibility.
The tested site RF is also a candidate, but the broad claim that nonlinear
methods cannot help is not supported by this small model grid. We have not fitted
interactions in the GAM or combined a nonlinear trait model with PIP.

These are extant-site transfer results with the existing phylogeny; a species can
occur at other training sites. Fossils have coarser, dated placements, and the
appropriate correction may differ. The next transfer check should deliberately
coarsen test taxonomy within each training/validation split. The choice of model
family is exploratory after comparing these outer scores; independent or further
prespecified validation is needed before treating the winning performance as a
confirmed production improvement. No fossil predictions or production fits were
replaced, and the old uncertainty formulas have not been validated for the new
recalibrated estimator.

## Reproduce and inspect

```
Rscript code/03f_map_model_experiment.R
Rscript code/03g_map_experiment_results.R
```

Requires mgcv and ranger plus existing project dependencies. `DILP_SOURCE` can
load a local checkout; `PIP_EXPERIMENT_R_LIB` adds a dependency library. This run
installed ranger in `/tmp/pip-map-r-library` and used the local dilp checkout.
`PIP_MAP_FOLDS` accepts disjoint outer-fold batches. Each checkpoint records
session/package versions and all inner predictions and memberships.

- `tables/map_model_experiment/model_comparison.csv`: full metrics.
- `tables/map_model_experiment/site_predictions.csv`: paired held-out predictions.
- `tables/map_model_experiment/fold_metrics.csv`: per-fold RMSE.
- `tables/map_model_experiment/recalibration_coefficients.csv`: ten corrections.
- `tables/map_model_experiment/selected_settings.csv`: inner-selected GAM/RF settings.
- `models/map_model_experiment/outer_XX.rds`: audit checkpoints.
- `plots/map_model_experiment_comparison.png`: PIP, recalibration, best tested GAM,
  and best tested RF, sharing axes.
- `plots/map_model_experiment_all_models.png`: all nine primary candidates.

## Exploratory fossil application

`code/04c_map_recalibration_trial.R` fits the final correction to all 92 extant
out-of-fold PIP site predictions (the nested experiment above assessed the
procedure). The resulting correction is

```
corrected ln(MAP cm) = -7.035270 + 2.358848 * original ln(MAP cm)
```

It is applied to the unrounded fossil site log predictions, followed by
exponentiation to cm. It leaves a prediction of about 177.2 cm unchanged,
reduces predictions below that value, and increases predictions above it.
The final slope differs from the mean inner-fold slope because this final
regression uses all 92 extant out-of-fold predictions.

Formal-taxonomy results:

| Site | Current MAP (cm) | Recalibrated MAP (cm) | Change |
|---|---:|---:|---:|
| Fox Hills | 154.1 | 127.4 | -17.3% |
| Williston Basin I | 169.5 | 159.5 | -5.9% |
| Palacio de los Loros | 166.3 | 152.5 | -8.3% |
| Williston Basin II | 162.1 | 143.6 | -11.4% |
| Williston Basin III | 160.8 | 140.8 | -12.4% |
| Cerrejon | 216.4 | 283.9 | +31.2% |
| Hubble Bubble | 166.4 | 152.8 | -8.2% |
| Laguna del Hunco | 167.8 | 155.7 | -7.2% |
| Republic | 146.5 | 113.0 | -22.8% |
| Bonanza | 148.6 | 116.9 | -21.3% |

Both taxonomy scenarios are saved in `tables/fossil_map_recalibration_trial.csv`.
Coefficients are in `tables/fossil_map_recalibration_coefficients.csv`; the plot
is `plots/fossil_map_recalibration_trial.png`. This is a separate point-estimate
sensitivity analysis. MAT and the production prediction tables are unchanged.
No fossil interval was recalibrated: these would need to account for the new
estimator and its correction uncertainty. Accuracy of this correction under
fossil taxonomic placements and ages has not yet been validated.

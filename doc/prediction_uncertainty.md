# Uncertainty in fossil climate predictions

## September 22, 2026 trial: site intervals fail coverage

The joint conditional PIP uncertainty calculation is implemented, but its nominal 95% intervals are **not calibrated for actual site climate**. Ten-fold site-grouped validation, refitting all models and imputers within training folds, gives:

| Target | Sites | Covered | Coverage | Mean interval width | Point RMSE |
|---|---:|---:|---:|---:|---:|
| MAT (degrees C) | 92 | 35 | 38.0% | 3.801 | 3.452 |
| ln(MAP in cm) | 92 | 29 | 31.5% | 0.401 | 0.529 |

MAP coverage here uses all 92 PIP sites. The three-model comparison uses 90 common eligible MAP sites and has PIP RMSE 0.535. These are different evaluation sets, not conflicting estimates.

Individual occurrence intervals cover local site climate in 93.9% of represented-taxon MAT cases and 94.1% of unrepresented-taxon cases; the MAP figures are 91.2% and 93.4%. Occurrences within a site are not independent validation replicates. These diagnostics do not validate intervals for species' grand-mean climate associations. Site coverage remains poor in both species-count groups, both training-representation groups, both missing-trait groups, and both halves of the observed climate range. The calculation loses too much uncertainty when aggregating to site climate.

**Do not replace fossil error bars with these intervals.** The outputs are conditional model intervals, retained to diagnose the failure. Current RMSE remains a descriptive performance benchmark, not a fossil-specific interval.

## Implemented calculation

For regression design matrices X_e and X_f, calibration covariance V_ee, fossil covariance V_ff, and fossil-to-calibration covariance C, define K = inverse(V_ee) and A = X_f - C K X_e. The joint prediction-error covariance is

```
sigma2 = residuals' K residuals / (n_training - n_coefficients)
Var(beta_hat) = sigma2 * inverse(X_e' K X_e)
S = sigma2 * (V_ff - C K C') + A Var(beta_hat) A'
```

The second term accounts for coefficient estimation after the PIP residual correction. Using X_f Var(beta_hat) X_f' instead would omit that adjustment. The code checks that S is positive semidefinite and fails on a material inconsistency rather than clipping invalid eigenvalues.

An individual occurrence's standard error is sqrt(S_ii). For an equally weighted assemblage of n occurrences, the site standard error is sqrt(sum(S) / n^2), including all off-diagonal terms. Approximate 95% Student intervals use residual degrees of freedom and treat the estimated lambda as fixed. For MAP, prediction, propagation, and intervals use natural-log units; exponentiated endpoints accompany the geometric site estimate in cm.

The implementation conditions on the fitted covariance parameters, measured or imputed traits, taxonomic placements, and point ages. It does not quantify uncertainty in lambda, imputation, calibration sampling, ages, or placement selection. Validation demonstrates inadequate site-climate coverage, but does not identify a unique missing error component.

The framework follows the PIP covariance approach discussed by [Gardner et al. (2025)](https://www.nature.com/articles/s41467-025-61036-1). Its application to an average of species predictions as an estimate of site climate requires the empirical check above.

## Covariance root correction

The trial exposed an existing inconsistency: calibration used the VCV of a pruned tree with root depth 152.8683 Ma, while fossil cross-covariances came from a scaffold with root depth 423 Ma. Every entry of the extant covariance block differed by about 270.1317 Ma (apart from numerical diagonal jitter). Those blocks cannot form the joint model required for prediction variance.

Climate fitting, validation, and fossil prediction now subset covariance matrices from the same rooted scaffold. Complete-case fits also preserve that root. Models were refitted, and the invariant is checked both in the fossil script and in `tests/test_climate_covariance_root.R`. Lambda transformation still scales off-diagonals only; diagonal jitter is 1e-6. Fossil prediction does not prune away the scaffold root before calculating covariance.

For a taxon represented at a training site, held-out prediction shares its phylogenetic component with training but uses a separate occurrence-specific residual. This matches the existing point predictor's lambda-scaled cross-covariance, including its same-taxon entries. Held-out taxa absent from the calibration scaffold retain the existing regression-only point prediction; their uncertainty uses root-attached contemporary tips with zero cross-covariance. These cases are flagged in validation output. This is an explicit fallback assumption, not evidence that their placements are known.

The related degradation and adjustment-field scripts now use the same fitted covariance scale. The taxonomic-degradation script also now uses node depths from the scaffold root for covariance, rather than ages before present returned by `branching.times()`.

## Files and reproduction

Run climate fitting (`02_phy_regression.R`), validation (`03_loso_cv.R`), coverage diagnostics (`03d_uncertainty_diagnostics.R`), June fossil cleaning (`00c_fossil_data_cleaning.R`), and fossil prediction (`04_fossil_predictions.R`) with the usual prerequisites.

- `tables/pip_cv_site_uncertainty.csv`: held-out site means, SEs, intervals, observations, coverage, and assemblage diagnostics.
- `tables/pip_cv_species_uncertainty.csv`: individual occurrence intervals, scaffold status, and training representation.
- `tables/pip_cv_interval_coverage.csv`: overall coverage and exploratory median-split strata.
- `tables/pip_cv_species_coverage.csv`: occurrence intervals compared with local site climate.
- `plots/pip_cv_conditional_intervals.png`: held-out observations and intervals.
- `tables/fossil_species_uncertainty_<scenario>.csv`: fossil occurrence intervals on MAT/log-MAP scales and response-scale endpoints.
- `tables/pip_cv_calibrated_coverage.csv`, `tables/pip_cv_calibrated_site_intervals.csv`: out-of-fold coverage of the calibrated intervals.
- `tables/fossil_site_uncertainty_<scenario>.csv`: conditional site intervals (`se`, `lower`, `upper`), a comparison SE that incorrectly assumes independent species errors, and the calibrated intervals (`site_discrepancy_sd`, `total_se`, `calibrated_lower`, `calibrated_upper`, and `_response` endpoints in cm for MAP). Report the calibrated interval.
- `models/fossil_prediction_covariance_<scenario>.rds`: named joint matrices for MAT and log(MAP).

Scenarios are `formal_only` and `include_informal`. The latter is a sensitivity analysis, not a probability distribution over placements. All new interval tables are separate from the existing fossil comparison tables. June source data are retained unchanged; the renamed comment header is normalized only during import.

Tests cover the independent-response limit, the zero-error identical-response limit, Monte Carlo covariance under repeated GLS fitting, row-order invariance, invalid covariance rejection, root consistency, and agreement of saved site intervals with full joint covariance. Existing placement, taxonomy, site-grouping, and site-local prediction checks also pass.

## Calibrated site intervals (shared site discrepancy)

The conditional intervals fail because most held-out site error is shared by every species at the site. In the full-data fit, residuals of single-site species at the same site correlate at about 0.85 (MAT) and 0.89 (log MAP); the phylogenetic covariance implies about 0.55 and 0.45, and after conditioning on training species almost none of that remains in the prediction-error covariance. Averaging species therefore divides away error that does not average out.

The calibrated interval adds a shared site discrepancy variance tau^2 to the conditional site variance:

    total_se = sqrt(tau^2 + se^2),   interval = estimate +/- qnorm(0.975) * total_se

tau^2 is estimated by moments from the 92 out-of-fold site residuals of `03_loso_cv.R`: tau^2 = mean(r^2) - mean(se^2), floored at zero. Mean bias (predictions slightly high: 0.62 degrees C for MAT, 0.16 log units for MAP) is absorbed into tau^2 rather than corrected. Current values: tau = 3.29 degrees C (MAT) and 0.518 (log MAP).

`03d_uncertainty_diagnostics.R` checks coverage with tau^2 estimated on the other nine folds only, so the reported coverage is out of sample: 96.7% for MAT and 93.5% for log MAP (conditional intervals: 38.0% and 31.5%). MAP errors have slightly heavier tails than normal.

For fossils, the conditional SE is much larger than for extant sites (2.2 to 4.2 degrees C, versus about 1 degree C) because deep, dated placements share little covariance with the calibration taxa. It is kept in the total rather than replaced by the RMSE. In CV the conditional SE carried little information about which sites erred (rank correlation with |residual| about 0.2), so how much it should widen fossil intervals is not validated. Coverage applies to sites exchangeable with the extant calibration sites; fossils add deep time, extinct lineages, and coarse placement, so their calibrated intervals are best treated as a lower bound on true uncertainty.

Before claiming fossil calibration, also test deliberately coarsened taxonomic placements and quantify sensitivity to training-site resampling and plausible imputations. Age uncertainty requires supplied age ranges or distributions. These extensions have not been run in this first trial.

## Independent coverage audit

`tests/audit_cv_coverage.R` independently reconstructed all 184 intervals from raw held-out traits and saved fold fits using direct linear-predictor error variance. It reproduces 35/92 MAT and 29/92 log(MAP) coverage, with SE agreement within 2e-11. See `doc/audits/cv_coverage_audit.md`. This checks the implementation, not whether the chosen covariance assumptions adequately describe local site-climate errors.

## Leave-one-site-out check (24 September 2026)

Refitting on 91 sites and predicting the remaining site, repeated across all 92,
gives conditional 95% coverage of **35/92 (38.0%) for MAT** and **25/92 (27.2%)
for log(MAP)**. Ten-fold values were 35/92 and 29/92. Point RMSE and mean interval
widths barely change. Leaving out fewer sites therefore does not resolve this
undercoverage. This evaluates the original conditional intervals, without a
site-discrepancy adjustment. See [the full check](audits/leave_one_site_out_coverage.md)
for validation details and reproduction commands.

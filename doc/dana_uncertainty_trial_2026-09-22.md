# Dana uncertainty trial — September 22, 2026

The conditional species-to-site uncertainty calculation is implemented and tested. It **fails the held-out site-climate coverage test** and should not supply published fossil error bars.

| Nominal 95% interval | Covered sites | Coverage | Average width |
|---|---:|---:|---:|
| MAT | 35 / 92 | 38.0% | 3.80 degrees C |
| ln(MAP in cm) | 29 / 92 | 31.5% | 0.401 log units |

Individual species intervals cover local climate roughly 91–94% of the time, but averaging their joint prediction errors gives inadequate site intervals. Phylogenetic and coefficient covariance alone are insufficient for site climate. The plot also shows systematic shrinkage toward intermediate climates. This result does not establish the exact missing error model.

## Completed

- Switched the fossil input to Dana's June file, preserving the original bytes and normalizing the comment header during import. The calibration file provided in the project root matches the active input exactly. There are still 360 fossil species-site occurrences across ten sites.
- Preserved the existing 92-site grouping and revised fossil ages.
- Implemented full prediction-error covariance, approximate conditional 95% intervals, arithmetic MAT averaging, and geometric MAP averaging. Intervals and covariance matrices are saved separately from existing comparison tables.
- Found and fixed inconsistent covariance roots: calibration previously used a 152.8683 Ma pruned tree while fossil prediction used a 423 Ma scaffold. All climate covariance blocks now preserve the scaffold root. Complete-case fits and related degradation/adjustment diagnostics use the same scale. The taxonomic-degradation diagnostic now uses node depths rather than ages before present for covariance.
- Refit climate PGLS models, reran all ten validation folds, regenerated both fossil taxonomy scenarios, refreshed climate validation/degradation/adjustment figures, and regenerated staged package model objects under `dilp_update/`.
- Synchronized affected release code, source data, and tests.

Current MAT PIP RMSE is 3.45154 degrees C. MAP PIP RMSE is 0.52949 across all 92 sites, or 0.53529 on the 90 sites eligible for the common DiLP comparison. These point-performance results are essentially unchanged; the failure is in interval coverage.

## Exploratory fossil MAT intervals

These intervals are conditional model outputs, **not validated site-climate intervals**. They condition on the supplied age, selected placement, completed traits, and fitted lambda. They exclude uncertainty in those inputs and parameters. Formal-only results:

| Site | Species | MAT estimate | Conditional lower | Conditional upper |
|---|---:|---:|---:|---:|
| Fox Hills | 24 | 16.46 | 10.67 | 22.24 |
| Williston Basin I | 20 | 17.30 | 9.35 | 25.24 |
| Palacio de los Loros | 29 | 17.17 | 12.77 | 21.58 |
| Williston Basin II | 23 | 16.52 | 10.06 | 22.98 |
| Williston Basin III | 18 | 16.49 | 11.30 | 21.67 |
| Cerrejon | 45 | 21.61 | 15.25 | 27.97 |
| Hubble Bubble | 16 | 19.16 | 14.07 | 24.26 |
| Laguna del Hunco | 119 | 16.66 | 8.49 | 24.83 |
| Republic | 41 | 14.09 | 9.03 | 19.15 |
| Bonanza | 25 | 17.10 | 11.57 | 22.64 |

Within-site covariance matters substantially for fossils. For example, Laguna del Hunco's conditional MAT SE is 4.16 degrees C; assuming independent species would give only 0.75 degrees C. This confirms that species count alone is not an adequate uncertainty measure, but does not establish interval calibration.

## Verification and limits

Analytical-limit and Monte Carlo GLS checks pass, as do covariance-root consistency, real fossil placement, site-local prediction, taxonomy, site grouping, aggregation, and generated-output checks. Both taxonomy scenarios retain all 360 occurrences. The coverage figure was inspected. Existing input outlier/margin warnings remain; the old PGLS summary routine also emits numerical F-test warnings. The prediction-covariance calculation and saved interval checks complete successfully.

The source checkout currently available at `/Users/jboyko/dilp` exposes the older `dilp_pgls(specimen_data)` API. A package parity attempt using a temporary copy with the refreshed model objects fails because that function does not accept `taxonomy_scenario`. The package implementation was not updated, installed, or published. Refreshed `dilp_update/` objects are staged artifacts, not a validated package release. Existing talk/manuscript Word files and the LMA analysis were not refreshed by this climate uncertainty trial.

The explicit root-placement fallback for extant taxa missing from the calibration scaffold is recorded in the validation tables. No separate fossil-like taxonomic-coarsening experiment, training-site bootstrap, repeated imputation, or age-range analysis was run. The two taxonomy scenarios are sensitivity cases, not probabilistically weighted uncertainties.

## Next test

Fit a shared site-discrepancy component using inner validation inside each outer training fold, then assess coverage on untouched outer sites. A common site-error variance would remain in the site interval rather than disappearing as more species are averaged. Assess systematic bias alongside interval width. Do not tune an inflation factor on these 92 held-out errors and present its coverage on the same errors as independent validation.

The implementation, equations, output definitions, and reproduction commands are in `doc/prediction_uncertainty.md`. The diagnostic figure is `plots/pip_cv_conditional_intervals.png`; numerical coverage is in `tables/pip_cv_interval_coverage.csv`.

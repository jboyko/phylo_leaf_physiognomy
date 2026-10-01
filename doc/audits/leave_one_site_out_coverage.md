# Leave-one-site-out coverage check — 24 September 2026

All 92 extant sites were held out individually. For each fit, raw specimens from
the other 91 sites were aggregated into training species means, training-side
imputers and regression/covariance parameters were refitted, and local held-out
species predictions were averaged to predict the excluded site's climate.
The imputed PIP model and conditional interval formula match the ten-fold trial.
The scaffold phylogeny is retained; species may occur at other training sites.

Coverage means the observed site MAT or log(MAP) lies between the nominal 95%
interval endpoints. Each site contributes one coverage indicator per target.

| Scheme | Target | Covered / sites | Coverage | RMSE | Mean interval width |
|---|---|---:|---:|---:|---:|
| ten_fold | mat | 35 / 92 | 38.0% | 3.4515 | 3.8012 |
| ten_fold | log_map | 29 / 92 | 31.5% | 0.5295 | 0.4012 |
| leave_one_site_out | mat | 35 / 92 | 38.0% | 3.4623 | 3.7835 |
| leave_one_site_out | log_map | 25 / 92 | 27.2% | 0.5327 | 0.4000 |

Leaving out one site instead of roughly nine does not resolve undercoverage.
MAT coverage is unchanged, MAP coverage is slightly lower, and point RMSE and
mean interval widths barely change. This comparison does not identify the cause
of undercoverage. It evaluates the original conditional intervals, not intervals
augmented with a fitted site-discrepancy variance or an RMSE-based error bar.

## Verification and reproduction

`code/03e_loso_coverage_comparison.R` verifies all 184 site-target rows, exactly
one held-out site and 91 training sites per fit, absence of overlap, observed
climate against `data/dat_site.csv`, site predictions against the mean of saved
species predictions, and interval endpoints against saved SEs and the training
residual degrees of freedom. Coverage flags are recomputed from endpoints.
All checks passed; all 92 fits completed without fitting/prediction errors.
The script also refuses to compile an incomplete set of checkpoints.

```sh
PIP_CV_SCHEME=leave_one_site_out Rscript code/03_loso_cv.R
Rscript code/03e_loso_coverage_comparison.R
```

This execution saved folds 1–8 serially and ran remaining folds in four disjoint
R-process batches via `PIP_CV_FOLDS`. Every fold uses seed `42 + fold`, as in the
serial procedure. Compact checkpoints preserve training/held-out site IDs,
training species, beta, lambda, and species/site predictions and intervals.
They omit dense covariance matrices; the earlier independent covariance audit
is documented separately in `cv_coverage_audit.md`.

Outputs are under `tables/leave_one_site_out/`: `coverage_comparison.csv`,
`paired_site_comparison.csv`, `validation_manifest.csv`, and the combined
`pip_cv_site_uncertainty.csv` and `pip_cv_species_uncertainty.csv`.
Plot: `plots/pip_leave_one_site_out_coverage.png`.
The default ten-fold outputs were preserved.

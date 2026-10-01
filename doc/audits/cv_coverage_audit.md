# Independent audit of reported PIP interval coverage

The reported counts are reproduced: **35 of 92 MAT site intervals** and **29 of 92 log(MAP) site intervals** contain the observed site climate. No production model or prediction outputs were changed during this audit.

## Checks

1. Independently counted interval inclusion in Python from numeric endpoints, without using saved coverage flags. Matched each record by site and target to the point-prediction table and the observed-climate input. Verified unique site-target records, fold labels, and the 92-site denominator for each target. Exponentiating log(MAP) endpoints and observations gives exactly the same inclusion decisions.
2. Reconstructed the intervals for all ten folds in R, using raw held-out morphotype traits and training-only imputers refitted with the recorded seeds. The audit reuses the production data-preparation definitions, but does not call the production uncertainty-fit, covariance, or interval functions.
3. Matched observed MAT and MAP against independent aggregation of the processed calibration source, rather than relying only on the observations copied into the interval table.
4. Reconstructed calibration covariance from the rooted scaffold and fitted lambda, and independently inverted it through Cholesky factorization. Checked against the fitted model's covariance and stored inverse. Reconstructed GLS coefficients and residual scale and checked the residual scale against the fitted PGLS RMS.
5. Expressed each species prediction as H times the training response vector Y. For the equally weighted site prediction, h is the mean of the rows of H. Computed its prediction-error variance directly as:

```
sigma2 * [mean(V_fossil) + h' V_training h - 2 h' mean_rows(C)]
```

This expands the variance of the difference between the future response mean and its linear predictor. It does not reuse the conditional-covariance-plus-coefficient-variance implementation. Checked H X_training = X_new, then independently reconstructed predictions, SEs, Student interval endpoints, and inclusion decisions.

## Agreement

- All 184 site-target predictions and intervals match.
- Maximum absolute SE difference: approximately 1.92e-11.
- Maximum absolute point-prediction difference: approximately 1.83e-11.
- MAT: 35 / 92 = 38.04%.
- log(MAP): 29 / 92 = 31.52%.

Reproduction: run `tests/audit_cv_coverage.R` from the project root, optionally setting `DILP_SOURCE` to a local dilp checkout. Numerical results are saved in `doc/audits/cv_coverage_independent.csv`. The audit requires the existing fold model files and calibration inputs.

## What this does and does not establish

This audit finds no counting, row-alignment, scale-conversion, variance-to-SE, or extra sample-size-division error in the implemented procedure. The independent reconstruction is conditional on the same chosen covariance model, including the separate occurrence residual for taxa represented in training and the root fallback for taxa absent from the calibration scaffold.

It does **not** establish that these covariance assumptions are the appropriate ones for local site climate. Nor does low coverage identify a unique cause. Bias, variance-model misspecification, omitted input uncertainty, the species-mean/site-climate distinction, or another mismatch could contribute. A shared site-error component remains a hypothesis, not a finding established by these counts. The result concerns this implementation and estimand; it is not proof that phylogeny-based uncertainty propagation generally fails.

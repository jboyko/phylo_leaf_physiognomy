# Uncertainty in fossil climate predictions

Cross-validation RMSE measures the typical prediction error across held-out extant sites. Each site is predicted by a model fitted without that site's measurements. The site errors are squared, averaged, and square-rooted:

\[
\mathrm{RMSE}=\sqrt{\frac{1}{N}\sum_{s=1}^{N}(\widehat y_s-y_s)^2}.
\]

The current imputed PIP MAT RMSE is 3.451 °C across 92 held-out site predictions. It measures predictive error, not the standard error of a fitted coefficient or a fossil site's estimated climate.

## Species and site uncertainty

PIP can estimate a predictive variance for each fossil occurrence, conditional on its traits, phylogenetic placement, age, and fitted model. The framework relates predictive uncertainty to phylogenetic covariance and separation from calibration taxa; see [Gardner et al., Phylogenetically informed predictions](https://www.nature.com/articles/s41467-025-61036-1).

For n species with equally weighted predictions, propagate their modeled prediction errors as

\[
\mathrm{Var}\left(\frac{1}{n}\sum_i \epsilon_i\right)
=\frac{1}{n^2}\sum_i\sum_j S_{ij},
\]

where epsilon_i is species i's prediction error and S is the joint prediction-error covariance matrix under the fitted model. Its diagonal contains species prediction-error variances. Its off-diagonal entries account for shared uncertainty, including phylogenetic relationships and estimated coefficients. More species can improve precision, but closely correlated errors reduce that gain. This quantifies uncertainty in the mean modeled species response; its calibration for the site's actual climate must be checked on held-out sites.

An analytical starting point, conditional on the fitted covariance parameters and measured or imputed traits, is

\[
S=\sigma^2(V_{ff}-CV_{ee}^{-1}C^T)
 + A\,\mathrm{Var}(\widehat\beta)\,A^T,
\qquad A=X_f-CV_{ee}^{-1}X_e.
\]

Here V_ee and V_ff are the lambda-transformed extant and fossil covariance blocks; C is the fossil–extant block; X_e and X_f are their regression design matrices. The residual scale is sigma squared. This is the prediction-error covariance for the modeled fossil responses; treating their mean as site climate requires empirical validation. The current climate output script calculates point predictions and does not yet report this matrix or site-specific prediction intervals.

For MAP, compute uncertainty for the mean log prediction, then exponentiate the interval endpoints to match the geometric site estimate.

## Validation and additional uncertainty

The first test is interval coverage on held-out extant sites: for example, whether nominal 95% intervals contain the observed site climate approximately 95% of the time. Check coverage across site composition, missing-data levels, and climate range, as well as overall coverage. Species in a site share climate, and some occur in the training data; these features can make an unvalidated species-level variance model overconfident for site climate.

The analytical calculation conditions on a fitted model and a chosen placement. Further uncertainty can be assessed by repeating the complete prediction under training-site resampling, plausible trait imputations, and alternative taxonomic placements. Age uncertainty requires age ranges or distributions; the supplied point ages alone do not define those distributions. Repeated predictions must preserve covariance among all species within a site. Validate and report these conditional and sensitivity components explicitly rather than combining them as though they were independent.

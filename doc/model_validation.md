# PART 3: Model Validation

## Methods

### Cross-validation design

Predictive accuracy was evaluated using 10-fold cross-validation with sites as the held-out unit, stratified by MAT. The 92 extant training sites were ranked by MAT and assigned round-robin to folds, producing folds of approximately 9 sites each. In each fold, all sites in that fold were withheld in full and all models were refitted from scratch on the remaining sites. The held-out sites were then predicted using the refitted models and the held-out sites' leaf trait measurements. Every site is predicted exactly once across the 10 folds. A fossil site has no representation in the extant training data, and this design matches that condition by withholding each site entirely from the models used to predict it. Held-out sites can nevertheless contain extant species represented in other training sites, so this validates site transfer with known extant relationships rather than the additional uncertainty of fossil placement.

The LMA analysis uses the same 10-fold framework with sites as the held-out unit. The 108 extant sites with known LMA are ranked by site-mean log₁₀(LMA) and assigned round-robin to folds, producing folds of approximately 11 sites each. At each fold, the held-out sites' morphotype measurements are withheld in full and all models are refitted from scratch on morphotypes from the remaining training sites, aggregated to species-level grand means. Held-out site predictions use within-site species means for the held-out morphotypes, with predictions averaged on the log₁₀(LMA) scale. This design matches the fossil prediction setting, where each fossil site contributes site-specific trait measurements per species.

### Model configurations

Twelve model configurations were evaluated, representing six model types crossed with two missing-data strategies. DiLP regression uses the published three-predictor equation for each climate target, with coefficients refitted within the same site folds. Yasuni-ridgetop and Yasuni-upper slope are combined before aggregation and fold assignment.

Four of the six model types are linear models (LM) differing in how the training data are aggregated. LM species is fitted to extant species-level grand means: a practical calibration choice that averages cross-site trait variation within species. Held-out predictions retain within-site species means. This species-trained LM is currently evaluated in validation only; the retained fossil LM is fitted to site means. This loss of within-species climatic and morphological variation during fitting should be considered when interpreting species-level comparisons. At prediction time, each held-out species is predicted using the refitted coefficients and species predictions are aggregated to a site-level estimate. The remaining three LM variants are fitted directly to site-level means and differ only in their tooth-trait treatment. Both LM site (specimen) and LM site (sp+zero) average the `dilp` morphotype means within each site, so they are identical by construction in the current data. LM site (untoothed excluded) uses the same site aggregation but leaves tooth traits as missing for confirmed untoothed morphotypes so they are excluded from tooth-trait averages. This is a local aggregation configuration, distinct from applying published DiLP coefficients.

The current species-level calibration labels are operational labels. Most correspond to named species, but records without a species epithet are pooled under genus-level keys and can span multiple sites and morphotypes. These unidentified pools are not equivalent to verified named species and remain pending taxonomic review. The present validation and numerical tables use this established label convention; separating those pools would change the training data and require a new full cross-validation run.

The fifth model type is PGLS, a phylogenetic generalized least squares model fitted to extant species-level means with Pagel's $\lambda$ jointly estimated. Predictions for held-out sites are computed as $X\beta$ per species and averaged to the site using the same procedure as LM species.

The sixth model type is PIP, the full Phylogenetically-Informed Prediction model described in Part 1. At each fold, PGLS coefficients $\beta$ and the $\lambda$-transformed variance-covariance matrix $V_{lam}$ are estimated from the training species. For each held-out site, the cross-covariance matrix $V_{cross}$ is computed between the $n$ training species and the held-out site's $m$ species using the phylogenetic variance-covariance function applied to the rooted extant scaffold. Because all held-out species are extant, no grafting is required and every entry of $V_{cross}$ is read directly from the existing tree. The prediction $\hat{y} = X\beta + V_{cross}^T V_{inv} e$ is evaluated per held-out species and averaged to the site. All LM and PGLS models use the same 12 fossil-measurable predictors defined in Part 1.

Two missing-data strategies were applied. Under bag imputation (impute), the imputer is fitted within each cross-validation fold using only the training data, then applied to the training and held-out data for that fold. Under complete-case analysis (CC), model fitting uses complete training rows. PIP and PGLS still impute held-out predictors and therefore evaluate all 92 sites. LM species and site CC variants do not impute held-out inputs: incomplete species or sites are excluded, with the evaluated site counts reported below.

The univariate LMA comparison used for fossil prediction evaluates three model types — LM, PGLS, and PIP — fitted to species-level means using a single predictor, log₁₀(petiole metric), where petiole metric is PW²/A (petiole width squared divided by blade area). Dana Royer (pers. comm.) established that LMA scales linearly with PW²/A in log–log space. Because the model has only one predictor, bag imputation is not applicable and the univariate comparison uses complete cases. The PGLS and PIP components follow the same fitting procedure described above, with Pagel's λ estimated jointly with the regression coefficients.

### Evaluation metrics

Three quantities are reported for each model and target variable across the $n$ held-out sites.

Root mean squared error (RMSE) is the primary accuracy metric.

$$\text{RMSE} = \sqrt{\frac{1}{n} \sum_{i=1}^{n} (\hat{y}_i - y_i)^2}$$

where $\hat{y}_i$ is the cross-validation prediction for observation $i$ (site for all targets) and $y_i$ is the observed value. For climate models, target variables are MAT in degrees Celsius and log(MAP) in log centimetres; for the LMA model the target is the site-mean log₁₀(LMA) in log₁₀ g m⁻². For species-level MAP models, log predictions are averaged across species within a site before RMSE is calculated; exponentiating this average gives the geometric site MAP estimate. This evaluates the same site-level estimator used for fossil MAP predictions. Site-level LM predictions are already site-level and remain on the log(MAP) scale for evaluation.

Mean bias $\bar{b}$ is the average signed prediction error.

$$\bar{b} = \frac{1}{n} \sum_{i=1}^{n} (\hat{y}_i - y_i)$$

Positive values indicate systematic over-prediction. RMSE decomposes into mean squared bias and residual variance.

$$\text{RMSE}^2 = \bar{b}^2 + \sigma^2_r$$

where $\sigma^2_r = \text{RMSE}^2 - \bar{b}^2$ is the variance of prediction errors after removing the mean offset. This decomposition identifies whether prediction error is driven by a systematic directional offset or by site-to-site variability in the errors.

The regression slope of $\hat{y}$ on $y$ (ordinary least squares, with predicted as response and observed as predictor) measures the degree to which predictions span the observed range. A slope of 1 indicates predictions that scale one-for-one with observations. A slope below 1 indicates predictions that are compressed toward the center of the observed distribution.

## Results

### Predictive accuracy

All scores below use held-out site predictions. Site coverage varies with missing predictors. PGLS and PIP use the same fitted regression; PIP adds the phylogenetic residual correction. Comparing their errors isolates the contribution of that correction.

| Model | MAT RMSE (°C) | ln(MAP) RMSE | n sites MAT / MAP |
| --- | --- | --- | --- |
| PIP (impute) | 3.452 | 0.529 | 92 / 92 |
| PIP (CC) | 3.733 | 0.512 | 92 / 92 |
| LM site, untoothed excluded (impute) | 3.767 | 0.594 | 92 / 92 |
| LM site, untoothed excluded (CC) | 3.788 | 0.615 | 81 / 81 |
| Refitted DiLP | 3.891 | 0.575 | 92 / 90 |
| LM site, sp+zero (CC) | 4.257 | 0.704 | 83 / 83 |
| LM site, sp+zero (impute) | 4.269 | 0.688 | 92 / 92 |
| LM species (CC) | 5.233 | 0.650 | 83 / 83 |
| LM species (impute) | 5.260 | 0.627 | 92 / 92 |
| PGLS (impute) | 5.890 | 0.633 | 92 / 92 |
| PGLS (CC) | 6.081 | 0.645 | 92 / 92 |

The direct comparison of PIP, the 12-trait site regression, and DiLP uses the same eligible sites within each target. These scores are stored in `tables/dana_cv_comparison.csv`.

![Cross-validated prediction error across model configurations.](../plots/fig1_rmse_comparison.png)

### Paired validation comparisons

Each comparison uses the same held-out sites for both models. Differences are PIP minus baseline RMSE. The percentile intervals use 10,000 paired site bootstrap draws conditional on the saved predictions and fold assignments. They summarize variation across evaluated sites and omit uncertainty from refitting the models.

| Target | Baseline | n sites | Δ RMSE | 95% site-bootstrap interval |
| --- | --- | --- | --- | --- |
| MAT (°C) | Refitted DiLP | 92 | -0.440 | [-1.066, 0.192] |
| MAT (°C) | 12-trait site LM | 92 | -0.818 | [-1.602, -0.068] |
| MAT (°C) | Untoothed-excluded site LM | 92 | -0.315 | [-0.929, 0.282] |
| MAT (°C) | PGLS | 92 | -2.439 | [-2.776, -2.081] |
| ln(MAP) | Refitted DiLP | 90 | -0.040 | [-0.112, 0.026] |
| ln(MAP) | 12-trait site LM | 92 | -0.158 | [-0.377, 0.027] |
| ln(MAP) | Untoothed-excluded site LM | 92 | -0.064 | [-0.174, 0.035] |
| ln(MAP) | PGLS | 92 | -0.104 | [-0.123, -0.084] |
| log₁₀(LMA) | Univariate LM | 107 | -0.010 | [-0.017, -0.004] |
| log₁₀(LMA) | Univariate PGLS | 107 | -0.011 | [-0.017, -0.005] |

PIP had lower observed error across all listed comparisons. The climate contrasts with refitted DiLP and the untoothed-excluded site LM span zero. The LMA contrasts and the climate contrasts with PGLS favor PIP.

### Prediction diagnostics

A slope below one for predicted versus observed climate indicates a compressed prediction range. Mean bias is the average signed error. The residual variance uses divisor n, so RMSE² equals bias² plus residual variance.

#### MAT (°C)

| Model | n | Slope | Mean bias | Bias² | Residual variance | RMSE |
| --- | --- | --- | --- | --- | --- | --- |
| PIP (impute) | 92 | 0.585 | 0.621 | 0.385 | 11.528 | 3.452 |
| LM site, untoothed-excl (impute) | 92 | 0.763 | -0.069 | 0.005 | 14.185 | 3.767 |
| DiLP regression (CV) | 92 | 0.707 | 0.011 | 0.000 | 15.141 | 3.891 |
| LM site, sp+zero (impute) | 92 | 0.720 | 0.035 | 0.001 | 18.227 | 4.269 |
| LM species (impute) | 92 | 0.340 | 1.778 | 3.160 | 24.506 | 5.260 |
| PGLS (impute) | 92 | 0.230 | 1.952 | 3.809 | 30.887 | 5.890 |

#### log(MAP)

| Model | n | Slope | Mean bias | Bias² | Residual variance | RMSE |
| --- | --- | --- | --- | --- | --- | --- |
| PIP (impute) | 92 | 0.238 | 0.158 | 0.025 | 0.256 | 0.529 |
| DiLP regression (CV) | 90 | 0.233 | 0.000 | 0.000 | 0.331 | 0.575 |
| LM site, untoothed-excl (impute) | 92 | 0.293 | 0.025 | 0.001 | 0.352 | 0.594 |
| LM species (impute) | 92 | 0.103 | 0.217 | 0.047 | 0.346 | 0.627 |
| PGLS (impute) | 92 | 0.095 | 0.227 | 0.052 | 0.349 | 0.633 |
| LM site, sp+zero (impute) | 92 | 0.256 | 0.043 | 0.002 | 0.471 | 0.688 |

### LMA validation

The single-predictor models used for fossil reconstruction were evaluated on the same 107 sites.

| Model | n | Slope | Mean bias | Bias² | Residual variance | RMSE |
| --- | --- | --- | --- | --- | --- | --- |
| PIP | 107 | 0.693 | 0.035 | 0.001 | 0.017 | 0.136 |
| LM | 107 | 0.573 | 0.032 | 0.001 | 0.021 | 0.147 |
| PGLS | 107 | 0.592 | 0.032 | 0.001 | 0.021 | 0.147 |

### Fossil uncertainty

Cross-validation RMSE summarizes extant prediction error. Fossil climate intervals propagate the full joint occurrence covariance and add a shared site discrepancy variance estimated from extant residuals. See [the interval method](prediction_uncertainty.md).
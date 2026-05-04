# PART 3: Model Validation

## Methods

### Cross-validation design

Predictive accuracy was evaluated using 10-fold cross-validation with sites as the held-out unit, stratified by MAT. The 92 extant training sites were ranked by MAT and assigned to folds using a snake pattern, producing folds of approximately 9 sites each. In each fold, all sites in that fold were withheld in full and all models were refitted from scratch on the remaining sites. The held-out sites were then predicted using the refitted models and the held-out sites' leaf trait measurements. Every site is predicted exactly once across the 10 folds. A fossil site has no representation in the extant training data, and this design matches that condition by withholding each site entirely from the models used to predict it.

### Model configurations

Twelve model configurations were evaluated, representing six model types crossed with two missing-data strategies.

Four of the six model types are linear models (LM) differing in how the training data are aggregated. LM species is fitted to extant species-level means. At prediction time, the held-out site's leaf traits are averaged by species, each species is predicted using the refitted coefficients, and species predictions are averaged to a site-level estimate. The remaining three LM variants are fitted directly to extant site-level means and differ only in how those means are computed from specimens. LM site (specimen) takes the mean of all specimens at each site regardless of species identity. LM site (sp+zero) first averages specimens within each species at each site, then averages across species, with tooth traits set to zero for confirmed untoothed species before averaging. LM site (Peppe) uses the same two-stage species-then-site aggregation but leaves tooth traits as missing for untoothed species so they are excluded from site-level tooth-trait averages, following Peppe et al. (2011).

The fifth model type is PGLS, a phylogenetic generalized least squares model fitted to extant species-level means with Pagel's $\lambda$ jointly estimated. Predictions for held-out sites are computed as $X\beta$ per species and averaged to the site using the same procedure as LM species.

The sixth model type is PIP, the full Phylogenetically-Informed Prediction model described in Part 1. At each fold, PGLS coefficients $\beta$ and the $\lambda$-transformed variance-covariance matrix $V_{lam}$ are estimated from the training species. For each held-out site, the cross-covariance matrix $V_{cross}$ is computed between the $n$ training species and the held-out site's $m$ species using the phylogenetic variance-covariance function applied to the full pruned extant tree. Because all held-out species are extant, no grafting is required and every entry of $V_{cross}$ is read directly from the existing tree. The prediction $\hat{y} = X\beta + V_{cross}^T V_{inv} e$ is evaluated per held-out species and averaged to the site. All LM and PGLS models use the same 12 fossil-measurable predictors defined in Part 1.

Two missing-data strategies were applied. Under bag imputation (impute), missing predictor values are filled using bagged-tree imputation fitted to the full species-level dataset before cross-validation begins, and the imputation model is not refitted at each fold. Under complete-case analysis (CC), only species with no missing values across the 12 predictors are retained, and sites that become empty after this filter are dropped from that fold's evaluation, reducing the effective sample to 90 sites for some configurations.

### Evaluation metrics

Three quantities are reported for each model and target variable across the $n$ held-out sites.

Root mean squared error (RMSE) is the primary accuracy metric.

$$\text{RMSE} = \sqrt{\frac{1}{n} \sum_{i=1}^{n} (\hat{y}_i - y_i)^2}$$

where $\hat{y}_i$ is the cross-validation predicted climate for site $i$ and $y_i$ is the observed climate. Target variables are MAT in degrees Celsius and log(MAP) in log centimetres.

Mean bias $\bar{b}$ is the average signed prediction error.

$$\bar{b} = \frac{1}{n} \sum_{i=1}^{n} (\hat{y}_i - y_i)$$

Positive values indicate systematic over-prediction. RMSE decomposes into mean squared bias and residual variance.

$$\text{RMSE}^2 = \bar{b}^2 + \sigma^2_r$$

where $\sigma^2_r = \text{RMSE}^2 - \bar{b}^2$ is the variance of prediction errors after removing the mean offset. This decomposition identifies whether prediction error is driven by a systematic directional offset or by site-to-site variability in the errors.

The regression slope of $\hat{y}$ on $y$ (ordinary least squares, with predicted as response and observed as predictor) measures the degree to which predictions span the observed range. A slope of 1 indicates predictions that scale one-for-one with observations. A slope below 1 indicates predictions that are compressed toward the center of the observed distribution.

## Results

### Predictive accuracy

PIP had the lowest RMSE for both MAT and log(MAP) under both missing-data strategies (Fig. 1). Among non-phylogenetic models, the site-level LM using the Peppe et al. (2011) aggregation performed best. For MAT, the species-level LM and PGLS had higher RMSE than all site-level LM variants. For log(MAP), LM species and PGLS outperformed the LM site (specimen) and LM site (sp+zero) configurations, though both remained worse than LM site (Peppe).

![Figure 1. Cross-validation RMSE for all twelve model configurations, faceted by target variable. Models are ordered by mean RMSE across both targets. Colour indicates model family.](../plots/fig1_rmse_comparison.png)

**Table 1.** Ten-fold cross-validation RMSE for all twelve model configurations, ordered by MAT RMSE. Values are root mean squared error for MAT (°C) and log(MAP) (log cm). Impute configurations have 93 held-out sites. CC configurations have 82 sites (LM site, Peppe) or 84 sites (all others) due to complete-case filtering.

| Model | MAT RMSE (°C) | log(MAP) RMSE | n sites |
| --- | --- | --- | --- |
| PIP (impute) | 3.417 | 0.525 | 93 |
| PIP (CC) | 3.709 | 0.508 | 93 |
| LM site, Peppe (CC) | 3.725 | 0.601 | 82 |
| LM site, Peppe (impute) | 3.736 | 0.590 | 93 |
| LM site, sp+zero (CC) | 4.207 | 0.683 | 84 |
| LM site, specimen (CC) | 4.207 | 0.683 | 84 |
| LM site, sp+zero (impute) | 4.209 | 0.672 | 93 |
| LM site, specimen (impute) | 4.217 | 0.675 | 93 |
| LM species (CC) | 5.186 | 0.644 | 84 |
| LM species (impute) | 5.198 | 0.621 | 93 |
| PGLS (impute) | 5.857 | 0.628 | 93 |
| PGLS (CC) | 6.054 | 0.640 | 93 |

PGLS and PIP share the same $\beta$ coefficients and phylogenetic covariance structure. PGLS predicts held-out sites from $X\beta$ alone. PIP adds the adjustment term $V_{cross}^T V_{inv} e$. PGLS has higher RMSE than the non-phylogenetic site-level LM for both targets, while PIP has lower RMSE. Comparing PGLS and PIP therefore isolates the predictive contribution of the phylogenetic adjustment term from the PGLS regression itself.

### Bias and residual variance

PIP had the lowest regression slope for both targets despite having the lowest RMSE (Fig. 2). Regression slopes and the bias-variance decomposition for all six impute models are given below.

![Figure 2. Observed versus predicted climate for three model configurations (impute variant). Rows are model types; columns are target variables. The dashed line is the 1:1 reference. The solid line and shaded band show the OLS regression of predicted on observed with a 95% confidence interval. Axes are independent per panel.](../plots/fig2_obs_vs_pred.png)

**Table 2.** Regression slope and RMSE decomposition for impute model configurations, MAT target (°C). Slope is the OLS regression slope of predicted on observed (a slope of 1 indicates predictions spanning the full observed range). Mean bias is the mean signed prediction error (positive values indicate over-prediction). RMSE² decomposes as Bias² + Residual variance.

#### MAT

| Model | Slope | Mean bias (°C) | Bias² | Residual var. | RMSE |
| --- | --- | --- | --- | --- | --- |
| PIP | 0.59 | +0.60 | 0.36 | 11.31 | 3.417 |
| LM site, Peppe | 0.78 | +0.04 | 0.00 | 13.95 | 3.736 |
| LM site, sp+zero | 0.73 | +0.05 | 0.00 | 17.72 | 4.209 |
| LM site, specimen | 0.73 | +0.05 | 0.00 | 17.78 | 4.217 |
| LM species | 0.35 | +1.73 | 3.00 | 24.02 | 5.198 |
| PGLS | 0.23 | +1.89 | 3.58 | 30.73 | 5.857 |

**Table 3.** Regression slope and RMSE decomposition for impute model configurations, log(MAP) target (log cm). Column definitions as in Table 2.

#### log(MAP)

| Model | Slope | Mean bias | Bias² | Residual var. | RMSE |
| --- | --- | --- | --- | --- | --- |
| PIP | 0.25 | +0.15 | 0.023 | 0.253 | 0.525 |
| LM site, Peppe | 0.31 | +0.02 | 0.000 | 0.347 | 0.590 |
| LM species | 0.11 | +0.21 | 0.044 | 0.342 | 0.621 |
| PGLS | 0.10 | +0.22 | 0.048 | 0.346 | 0.628 |
| LM site, sp+zero | 0.28 | +0.04 | 0.002 | 0.449 | 0.672 |
| LM site, specimen | 0.27 | +0.04 | 0.001 | 0.454 | 0.675 |

PIP has the lowest RMSE and the lowest regression slope for both targets. For MAT, PIP's slope is 0.59 versus 0.78 for the LM site (Peppe) model. For log(MAP), PIP's slope is 0.25 versus 0.31. PIP's lower slope indicates that its predictions do not span the observed climate range. Relative to the LM site (Peppe) model, PIP reduces residual variance by 2.64 °C² for MAT and by 0.094 for log(MAP), while increasing squared bias by 0.36 °C² for MAT and by 0.023 for log(MAP). The net change in RMSE² favors PIP in both cases.

The phylogenetic adjustment $V_{cross}^T V_{inv} e$ weights each training species' decorrelated residual by its phylogenetic proximity to the held-out species. Training residuals average near zero across the full dataset, so for held-out species with no close phylogenetic neighbors carrying large, consistent residuals the adjustment is small and pulls the prediction toward the training mean. For sites at the extremes of the MAT and MAP distributions this produces a systematic compression of predictions toward the center of the observed range, which is visible as the sub-unity slopes in Figure 2. Fossil site predictions that fall near the boundaries of the training climate distribution will carry this shrinkage, and point estimates for such sites should be interpreted accordingly.

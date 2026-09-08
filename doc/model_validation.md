# PART 3: Model Validation

## Methods

### Cross-validation design

Predictive accuracy was evaluated using 10-fold cross-validation with sites as the held-out unit, stratified by MAT. The 93 extant training sites were ranked by MAT and assigned round-robin to folds, producing folds of approximately 9 sites each. In each fold, all sites in that fold were withheld in full and all models were refitted from scratch on the remaining sites. The held-out sites were then predicted using the refitted models and the held-out sites' leaf trait measurements. Every site is predicted exactly once across the 10 folds. A fossil site has no representation in the extant training data, and this design matches that condition by withholding each site entirely from the models used to predict it. Held-out sites can nevertheless contain extant species represented in other training sites, so this validates site transfer with known extant relationships rather than the additional uncertainty of fossil placement.

The LMA analysis uses the same 10-fold framework with sites as the held-out unit. The 108 extant sites with both petiole metric measurements and known LMA are ranked by site-mean log₁₀(LMA) and assigned round-robin to folds, producing folds of approximately 11 sites each. At each fold, the held-out sites' morphotype measurements are withheld in full and all models are refitted from scratch on morphotypes from the remaining training sites, aggregated to species-level grand means. Held-out site predictions use within-site species means for the held-out morphotypes, averaged to a site-level estimate. This design matches the fossil prediction setting, where each fossil site contributes site-specific trait measurements per species.

### Model configurations

Twelve model configurations were evaluated, representing six model types crossed with two missing-data strategies.

Four of the six model types are linear models (LM) differing in how the training data are aggregated. LM species is fitted to extant species-level grand means: a practical calibration choice that averages cross-site trait variation within species. Held-out and fossil predictions retain within-site species means. This loss of within-species climatic and morphological variation during fitting should be considered when interpreting species-level comparisons. At prediction time, each held-out species is predicted using the refitted coefficients and species predictions are aggregated to a site-level estimate. The remaining three LM variants are fitted directly to site-level means and differ only in their tooth-trait treatment. Both LM site (specimen) and LM site (sp+zero) average the `dilp` morphotype means within each site, so they are identical by construction in the current data. LM site (untoothed excluded) uses the same site aggregation but leaves tooth traits as missing for confirmed untoothed morphotypes so they are excluded from tooth-trait averages. This is a local aggregation configuration, distinct from applying published DiLP coefficients.

The current species-level calibration labels are operational labels. Most correspond to named species, but records without a species epithet are pooled under genus-level keys and can span multiple sites and morphotypes. These unidentified pools are not equivalent to verified named species and remain pending taxonomic review. The present validation and numerical tables use this established label convention; separating those pools would change the training data and require a new full cross-validation run.

The fifth model type is PGLS, a phylogenetic generalized least squares model fitted to extant species-level means with Pagel's $\lambda$ jointly estimated. Predictions for held-out sites are computed as $X\beta$ per species and averaged to the site using the same procedure as LM species.

The sixth model type is PIP, the full Phylogenetically-Informed Prediction model described in Part 1. At each fold, PGLS coefficients $\beta$ and the $\lambda$-transformed variance-covariance matrix $V_{lam}$ are estimated from the training species. For each held-out site, the cross-covariance matrix $V_{cross}$ is computed between the $n$ training species and the held-out site's $m$ species using the phylogenetic variance-covariance function applied to the full pruned extant tree. Because all held-out species are extant, no grafting is required and every entry of $V_{cross}$ is read directly from the existing tree. The prediction $\hat{y} = X\beta + V_{cross}^T V_{inv} e$ is evaluated per held-out species and averaged to the site. All LM and PGLS models use the same 12 fossil-measurable predictors defined in Part 1.

Two missing-data strategies were applied. Under bag imputation (impute), the imputer is fitted within each cross-validation fold using only the training data, then applied to the training and held-out data for that fold. Under complete-case analysis (CC), model fitting uses complete training rows. PIP and PGLS still impute held-out predictors and therefore evaluate all 93 sites. LM species and site CC variants do not impute held-out inputs: incomplete species or sites are excluded, leaving 84 evaluated sites for LM species and zero-filled site configurations and 82 for the untoothed-excluded site configuration.

The LMA analysis evaluates three model types — LM, PGLS, and PIP — fitted to species-level means using a single predictor, log₁₀(petiole metric), where petiole metric is PW²/A (petiole width squared divided by blade area). Dana Royer (pers. comm.) established that LMA scales linearly with PW²/A in log–log space. Because the model has only one predictor, bag imputation is not applicable and a single complete-case configuration is used throughout. The PGLS and PIP components follow the same fitting procedure described above, with Pagel's λ estimated jointly with the regression coefficients.

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

PIP had the lowest cross-validated RMSE for both MAT and log(MAP) under both missing-data strategies. The headline comparison distinguishes PIP impute, the practical LM site sp+zero impute baseline, and fixed published DiLP coefficients. The latter is an in-sample reference, not a newly cross-validated competitor: its scores were not obtained by the same evaluation. PGLS versus PIP remains the strict ablation of the prediction-time phylogenetic correction; PIP versus LM site changes both the treatment of phylogeny and the fitting/aggregation strategy.

![Figure 1. Cross-validation RMSE for twelve model configurations plus the published in-sample reference, faceted by target variable. Models are ordered by mean RMSE across both targets. Colour indicates model family.](../plots/fig1_rmse_comparison.png)

**Table 1.** Ten-fold cross-validation RMSE for twelve model configurations, ordered by MAT RMSE, followed by the fixed published DiLP in-sample reference. Values are root mean squared error for MAT (°C) and log(MAP) (log cm). Impute configurations evaluate 93 held-out sites; PIP and PGLS complete-case configurations also evaluate 93 sites because held-out inputs are imputed.

| Model | MAT RMSE (°C) | ln(MAP) RMSE | n sites MAT / MAP |
| --- | ---: | ---: | ---: |
| PIP (impute) | 3.414 | 0.525 | 93 |
| PIP (CC) | 3.707 | 0.508 | 93 |
| LM site, untoothed excluded (impute) | 3.722 | 0.583 | 93 |
| LM site, untoothed excluded (CC) | 3.725 | 0.601 | 82 |
| LM site, sp+zero (CC) | 4.207 | 0.683 | 84 |
| LM site, specimen (CC) | 4.207 | 0.683 | 84 |
| LM site, sp+zero (impute) | 4.217 | 0.672 | 93 |
| LM site, specimen (impute) | 4.217 | 0.672 | 93 |
| LM species (CC) | 5.186 | 0.644 | 84 |
| LM species (impute) | 5.197 | 0.621 | 93 |
| PGLS (impute) | 5.854 | 0.628 | 93 |
| PGLS (CC) | 6.052 | 0.640 | 93 |
| Published DiLP (in-sample) | 3.791 | 0.554 | 93 / 91 |

PGLS and PIP share the same $\beta$ coefficients and phylogenetic covariance structure. PGLS predicts held-out sites from $X\beta$ alone. PIP adds the adjustment term $V_{cross}^T V_{inv} e$. PGLS has higher RMSE than the best-performing site-level LM for both targets, while PIP has lower RMSE. Comparing PGLS and PIP therefore isolates the predictive contribution of the phylogenetic adjustment term from the PGLS regression itself.

### Bias and residual variance

PIP had the lowest RMSE for both climate targets despite having compressed predictions, as measured by the regression slope of predicted on observed (Fig. 2). Regression slopes and the bias-variance decomposition for impute models are given below.

![Figure 2. Observed versus predicted climate for three cross-validated imputed configurations and the published in-sample reference. Rows are model types; columns are target variables. The dashed line is the 1:1 reference. The solid line and shaded band show the OLS regression of predicted on observed with a 95% confidence interval. Axes are independent per panel.](../plots/fig2_obs_vs_pred.png)

**Table 2.** Regression slope and RMSE decomposition for impute model configurations, MAT target (°C). Slope is the OLS regression slope of predicted on observed (a slope of 1 indicates predictions spanning the full observed range). Mean bias is the mean signed prediction error (positive values indicate over-prediction). RMSE² decomposes as Bias² + Residual variance.

#### MAT

| Model | Slope | Mean bias (°C) | Bias² | Residual var. | RMSE |
| --- | --- | --- | --- | --- | --- |
| PIP | 0.59 | +0.60 | 0.36 | 11.30 | 3.414 |
| LM site, untoothed excluded | 0.78 | +0.00 | 0.00 | 13.85 | 3.722 |
| LM site, sp+zero | 0.73 | +0.06 | 0.00 | 17.78 | 4.217 |
| LM site, specimen | 0.73 | +0.06 | 0.00 | 17.78 | 4.217 |
| LM species | 0.35 | +1.73 | 3.00 | 24.01 | 5.197 |
| PGLS | 0.23 | +1.89 | 3.56 | 30.70 | 5.854 |

**Table 3.** Regression slope and RMSE decomposition for impute model configurations, log(MAP) target (log cm). Column definitions as in Table 2.

#### log(MAP)

| Model | Slope | Mean bias | Bias² | Residual var. | RMSE |
| --- | --- | --- | --- | --- | --- |
| PIP | 0.25 | +0.15 | 0.023 | 0.253 | 0.525 |
| LM site, untoothed excluded | 0.31 | +0.02 | 0.000 | 0.339 | 0.583 |
| LM species | 0.11 | +0.21 | 0.044 | 0.341 | 0.621 |
| PGLS | 0.10 | +0.22 | 0.048 | 0.346 | 0.628 |
| LM site, sp+zero | 0.28 | +0.04 | 0.002 | 0.451 | 0.672 |
| LM site, specimen | 0.28 | +0.04 | 0.002 | 0.451 | 0.672 |

PIP's lower slope indicates that its predictions do not span the full observed climate range. The non-phylogenetic comparison is the LM site untoothed-excluded configuration; it must not be described as fixed published DiLP coefficients.

The phylogenetic adjustment weights training residuals by covariance with the held-out taxa. The sub-unity prediction slopes document compression of the reconstructed climate range, but do not establish that the correction always moves predictions toward the training mean. Fossil reconstructions near the limits of the calibration climate range should be interpreted with this observed limitation and the additional uncertainty in fossil placement in mind.

### LMA predictive accuracy

PIP had the lowest RMSE for log₁₀(LMA); PGLS and LM performed similarly to each other (Table 4). All three models were evaluated on 108 complete-case sites.

**Table 4.** Ten-fold cross-validation RMSE for LMA models. Values are root mean squared error for site-mean log₁₀(LMA) (log₁₀ g m⁻²). Folds are stratified by site-mean log₁₀(LMA).

| Model | log₁₀(LMA) RMSE | n sites |
| --- | --- | --- |
| PIP | 0.130 | 108 |
| PGLS | 0.144 | 108 |
| LM | 0.146 | 108 |

### LMA bias and residual variance

The bias-variance decomposition for site-mean log₁₀(LMA) is given in Table 5. All three models show a small positive mean bias of 0.02–0.03 log₁₀ g m⁻², and the differences across models are negligible. PIP's RMSE advantage over LM comes entirely from reduced residual variance (0.016 versus 0.020). This contrasts with the climate result, where PIP carried slightly higher bias than LM site (untoothed excluded) alongside lower variance. For LMA, there is no bias cost: the adjustment reduces variance without trading it against accuracy in the mean. PIP also has the highest regression slope (0.665) among the three LMA models, followed by LM (0.577) and PGLS (0.543), meaning PIP predictions span more of the observed LMA range than either alternative. For climate, PIP predictions were more compressed than site-level LM predictions, but less compressed than species-level LM and PGLS predictions.

**Table 5.** Regression slope and RMSE decomposition for LMA models, site-mean log₁₀(LMA) target. Column definitions as in Tables 2–3.

| Model | Slope | Mean bias | Bias² | Residual var. | RMSE |
| --- | --- | --- | --- | --- | --- |
| PIP | 0.665 | +0.028 | 0.001 | 0.016 | 0.130 |
| PGLS | 0.543 | +0.022 | 0.000 | 0.020 | 0.144 |
| LM | 0.577 | +0.032 | 0.001 | 0.020 | 0.146 |

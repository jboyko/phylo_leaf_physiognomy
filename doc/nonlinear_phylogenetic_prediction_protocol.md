# Protocol: nonlinear phylogenetic prediction of fossil MAP

## Objective and scope

Test whether a nonlinear relationship between fossil-measurable leaf traits and climate improves site-level MAP prediction while retaining prediction-time information from phylogenetic relatives.

Implement a kernel model with phylogenetically correlated residuals. This is an adaptation motivated by Rosas-Puchuri et al. (2024), **not a reproduction of their phyloKRR algorithm**, which applies a nonlinear kernel to phylogenetically transformed observations. Here the nonlinear function acts on the original, training-standardized traits. Do not whiten the complete response vector before splitting sites.

Start with log(MAP in cm). Preserve the existing production models, fossil estimates, and uncertainty tables. Save this experiment separately. The first deliverable is validated point prediction; fossil interval calibration is a separate question.

Background: [paper](https://doi.org/10.1111/2041-210X.14385), [reading notes](../summary.md), [previous model experiment](map_model_experiment.md).

## 1. Reuse the existing data rules

Use the current calibration input and preprocessing in `code/03_loso_cv.R` and the nested split pattern in `code/03f_map_model_experiment.R`.

1. Combine Yasuni-ridgetop and Yasuni-upper slope before aggregation and fold assignment. Expect 92 sites.
2. Apply the documented untoothed-leaf filling rules before aggregation.
3. Use only the authoritative `fossil_traits` predictors. Do not add observed climate, coordinates, site identity, or non-fossil-measurable traits as predictors.
4. For each training split, reconstruct species-level trait and climate means using only that split's sites. Retain the current operational taxon labels and aggregation conventions for comparability.
5. Retain the current training response: log of the species' aggregated MAP in cm. Do not silently change this to mean log(MAP).
6. Estimate predictor eligibility, imputation, centering, scaling, and constant-column removal using training data only. Apply the fitted transformations unchanged to validation/test rows.
7. For a held-out site, aggregate its traits within species at that site. Never pool its traits with occurrences at other sites.
8. Average occurrence predictions on the log scale to obtain the site prediction. Exponentiating gives the existing geometric MAP estimator in cm. Do not introduce a lognormal mean correction or average exponentiated occurrence predictions.

For a clean new implementation, use a trait-only training imputer for both training and prediction. Audit this against the existing PIP path, whose training PGLS imputer also includes the training response. If this changes completed training traits, report the difference and include a matched-preprocessing PIP control. Do not attribute a preprocessing change to the kernel.

## 2. Define a coherent joint covariance model

Use a Gaussian-process/kernel representation with a fixed regression mean:

    y = H beta + u_trait + u_phylogeny + epsilon

Here y contains training species log(MAP), H is the chosen mean-design matrix, u_trait supplies a regularized trait function, and u_phylogeny + epsilon supplies phylogenetic residual covariance. The components are independent in this model.

Let B be the Brownian shared-branch covariance taken from the existing rooted scaffold. Divide **every** training, cross, and prediction block by the same fixed scaffold root height to improve numerical scaling. Do not renormalize separately by fold, species subset, or fossil age. Keep biological branch variance separate from numerical diagonal jitter.

Define:

    V_lambda = lambda_phy * B + (1 - lambda_phy) * diag(diag(B))
    A = V_lambda + eta * K_trait
    Cov(y) = sigma2 * A

- `lambda_phy` is Pagel-type residual phylogenetic dependence, between 0 and 1.
- `eta` is the trait-kernel variance relative to the residual covariance, at least zero. It controls flexibility/regularization together with the kernel shape.
- `sigma2` is an overall variance scale; it cancels from point-prediction weights.
- Off-diagonals of B are multiplied by lambda_phy; its diagonal is unchanged in V_lambda.
- The diagonal residual component follows the existing project's branch-depth convention. It is not a new model of leaf measurement error or site sampling error.

Use distinct names for lambda_phy, eta, numerical jitter, and any kernel bandwidth. Do not confuse the paper's ridge penalty called lambda with Pagel's lambda.

### Trait kernels

After training-only scaling, let z_i denote a p-dimensional trait vector.

    Linear: K(i,j) = z_i' z_j / p
    RBF:    K(i,j) = exp(-||z_i-z_j||^2 / (2 * ell^2))

The RBF kernel permits curved effects and interactions. Compute the reference distance for ell from training traits only, using the median positive pairwise Euclidean distance. Retain the same scaling and ell for all prediction blocks. Fail clearly if there are no usable predictors or positive distances.

### Candidate definitions

Evaluate these separately:

| Candidate | H | Trait component | Purpose |
|---|---|---|---|
| Existing PIP | Existing linear design | None | Established benchmark |
| Linear kernel + phylogeny | Intercept only | Linear kernel | Regularized linear control |
| RBF kernel + phylogeny | Intercept only | RBF kernel | Nonlinear counterpart to the linear-kernel control |
| Linear mean + RBF + phylogeny | Intercept plus traits | RBF kernel | Nonlinear extension retaining the existing linear mean |
| Nested recalibrated PIP | Existing PIP plus trained site correction | None | Best prior benchmark |

Do not put a linear kernel alongside an unpenalized identical linear mean and claim this isolates a useful additional linear effect: the components overlap. The eta=0 version of the last kernel candidate gives the PIP equations for matching H, covariance, imputation, and lambda_phy. This is a required numerical check, not a guarantee that different hyperparameter-selection procedures produce identical fitted models.

## 3. Fit and predict without using unknown responses

For fixed covariance hyperparameters, estimate the unpenalized mean by GLS:

    beta_hat = (H' A^-1 H)^-1 H' A^-1 y
    alpha = A^-1 (y - H beta_hat)

Use Cholesky solves rather than explicit matrix inverses. Detect rank deficiency in H and use a deterministic training-derived full-rank design consistently at prediction time. Record any removed columns. Do not mask rank deficiency with a large ridge penalty on the mean.

For held-out occurrences, form a cross-covariance block:

    C = lambda_phy * B_new,train + eta * K_new,train
    prediction_new = H_new beta_hat + C alpha

This formula includes both the nonlinear trait prediction and the prediction-time phylogenetic adjustment **once**. Do not append the existing PIP adjustment again; that would double-count information.

Save the separate contributions for diagnosis:

    linear_mean       = H_new beta_hat
    trait_adjustment  = eta * K_new,train * alpha
    phylo_adjustment  = lambda_phy * B_new,train * alpha

A held-out occurrence of a species represented elsewhere shares its phylogenetic component with training. It does not share the independent residual/nugget with the training row. Thus C contains no independent-residual term even for the same taxon. This matches the current occurrence-prediction convention.

Predict each occurrence first, then average its predicted log(MAP) within site. No held-out climate enters imputation, covariance fitting, prediction, or averaging.

## 4. Use nested site-grouped validation

Use the existing ten outer site folds so all methods predict exactly the same 92 sites as the previous experiment. The ten-fold PIP benchmark is RMSE 0.5295 and nested recalibrated PIP is 0.4220 in log units. The separate leave-one-site-out benchmark uses different fits and is not the paired comparator here.

For each outer fold:

1. Remove all specimens from its held-out sites before constructing training species means.
2. Partition the remaining sites into five reproducible inner folds. Reuse the previous experiment's saved inner memberships where available.
3. For each inner fold, rebuild preprocessing and all training covariance/kernel blocks using its training sites only.
4. Fit candidate hyperparameters and predict the inner held-out occurrences, then their sites.
5. Select hyperparameters by pooled inner **site-level** mean squared error in log(MAP). Weight sites equally, not by leaf/species count. Do not select by species-level error, in-sample likelihood, or apparent prediction spread.
6. Rebuild preprocessing on all outer training sites, fit the selected model, and predict the outer test sites once.
7. Save predictions, selected parameters, training/held-out IDs, inner losses, preprocessing metadata, diagnostic contributions, and seeds.

Initial fixed grid:

- lambda_phy: 0, 0.5, 0.9, 0.99.
- eta: 0.1, 1, 10.
- RBF ell: 0.5, 1, 2 times the training reference distance.
- Include eta=0 separately for the linear-mean candidate to allow selection of no nonlinear component.

The linear kernel does not have a bandwidth. These are exploratory ranges, not established optimal settings. Record boundary selections. If expansion is needed, define it using training-side diagnostics and rerun the full affected nested procedure; do not choose ranges after inspecting outer errors. Document any subsequent exploration as such.

Cache fold-local preprocessing, pairwise trait distances, and tree blocks. Parallelize independent outer folds with deterministic seeds. Benchmark one outer fold before dispatching the full run; avoid an uncontrolled large hyperparameter search.

Recalibrated PIP must retain its existing nested training-only correction. Do not fit a correction to all outer predictions and score it against the same outcomes. Initially evaluate the new kernels without a further recalibration; any later kernel recalibration requires its own properly nested fit.

## 5. Required verification

Before interpreting real-data results, check:

- **PIP limit:** eta=0 with the matching linear mean and covariance reproduces existing PIP predictions, including cross-covariance adjustments, to numerical tolerance.
- **No phylogenetic dependence:** lambda_phy=0 removes phylogenetic cross-covariance, leaving the chosen trait model and independent residuals.
- **Permutation invariance:** jointly reordering calibration rows and covariance blocks leaves predictions unchanged; reordering query occurrences only reorders predictions.
- **Batch invariance:** predicting a fossil alone or alongside unrelated query fossils yields the same point estimate under a fixed fitted model and placements.
- **Joint covariance:** training/new/cross blocks use one root scale and form a positive-semidefinite joint covariance before numerical jitter. No arbitrary eigenvalue clipping.
- **Response exclusion:** changing an outer held-out climate cannot change that fold's fitted objects or predictions. Inspect this for imputation, hyperparameter selection, and recalibration.
- **Membership:** no outer test site appears anywhere in its inner training/validation data. Every outer site receives one prediction per candidate.
- **Aggregation:** saved site log predictions equal the mean of saved occurrence log predictions; exponentiation gives the saved MAP in cm.
- **Missing traits:** all completed predictors and predictions are finite. Reuse documented missing-taxonomy/scaffold fallbacks consistently rather than silently dropping occurrences.
- **Simulated nonlinear signal:** a small synthetic dataset with a known curved trait function and phylogenetic residuals verifies the implementation can recover that relationship. Also test a linear case. These tests are implementation checks, not evidence for improvement on leaves.

Record numerical jitter and investigate materially negative conditional variances or failed factorizations. Scale tolerances relative to matrix magnitude.

## 6. Report the experiment

For each candidate, report paired outer predictions for all 92 sites and:

- RMSE and MAE on log(MAP).
- Mean signed error.
- Predicted-versus-observed correlation.
- Ratio of prediction SD to observed SD.
- Per-outer-fold RMSE and the selected hyperparameters.
- Mean, trait-kernel, and phylogenetic prediction contributions.

Plot observed versus predicted values on common axes, with a 1:1 line. Check errors across the observed climate range descriptively. Increased spread alone is not an improvement. Do not treat all leaves/species as independent site-validation replicates or claim that ten overlapping CV fits provide ten independent experiments.

Report both comparisons: whether the nonlinear model improves on its matched linear-kernel control, and whether it improves on nested recalibrated PIP. Selecting the best model family after inspecting these scores remains exploratory even though hyperparameter selection is nested.

Suggested new outputs:

    models/map_kernel_experiment/outer_XX.rds
    tables/map_kernel_experiment/site_predictions.csv
    tables/map_kernel_experiment/occurrence_predictions.csv
    tables/map_kernel_experiment/model_comparison.csv
    tables/map_kernel_experiment/selected_parameters.csv
    tables/map_kernel_experiment/validation_manifest.csv
    plots/map_kernel_experiment_comparison.png
    doc/map_kernel_experiment.md

## 7. Test transfer toward fossils

Only after the extant comparison is understood, evaluate coarsened test taxonomy (genus/family/order) using training-only placement evidence and the existing placement rules. Perform this within held-out site folds. This checks part of the fossil problem; it does not validate deep-time ages or extinct lineages.

For a fossil sensitivity application:

1. Select the candidate and fitting protocol explicitly from the documented experiment.
2. Tune it using site-grouped CV on all extant data, then fit it to all extant training species.
3. Apply its saved imputer and scaling to each fossil species-by-site occurrence's local traits.
4. Use the occurrence's own site age and genus → family → order → root placement fallback. Resolve placements from the original extant scaffold; never use an earlier grafted fossil as taxonomic evidence.
5. Construct B_new,train from those placements using the same scaffold root and normalization as calibration. Kernel cross-covariance uses the occurrence's own completed/scaled traits.
6. Predict occurrences with the formula in Section 3 and aggregate within site on the log scale.
7. Save both formal-only and informal-taxonomy sensitivity results, with original PIP and recalibrated PIP comparators. Do not pool fossil traits or ages across sites.
8. Preserve Dana's site grouping and revised ages, including Palacio de los Loros at 64.08 Ma.

Flag trait-space extrapolation and weak kernel similarity to calibration data. An RBF component tends to lose influence far from its training examples; a wider range of predictions in extant CV does not guarantee useful fossil extrapolation.

## 8. Keep uncertainty claims separate

A coherent joint kernel model permits conditional prediction covariance. For fixed hyperparameters, define:

    A_new = V_lambda,new + eta * K_new,new
    D = H_new - C A^-1 H
    S = sigma2 * [A_new - C A^-1 C' + D (H' A^-1 H)^-1 D']

Here A_new includes the independent new-occurrence residual if predicting an occurrence response. Distinguish that target from a latent trait-function value. If a conditional variance scale is needed, the fixed-shape GLS estimate is `(y-H beta_hat)' A^-1 (y-H beta_hat)/(n-rank(H))`; it does not account for selecting covariance hyperparameters.

For the mean of n occurrence predictions, propagate the full covariance as `sum(S)/n^2`. This is **not automatically uncertainty for actual site climate**. Our previous conditional PIP intervals failed that target even with correct covariance propagation. This extension does not remove the species-average/local-site distinction.

Treat kernel-based intervals as diagnostic until their held-out site coverage is measured. Do not reuse the old PIP Student degrees of freedom or attach its old RMSE to new predictions without evaluation. Hyperparameter fitting, imputation, specimen sampling, uncertain placements, and ages introduce uncertainty beyond this conditional formula. Any empirical interval calibration must also be fitted without access to its evaluation sites.

## Completion criteria

Deliver reproducible code, passing model-limit and leakage checks, paired nested predictions for all 92 sites, the full comparison and plots, and a written conclusion separating measured improvements from proposed explanations. A negative result is a completed experiment. Fossil sensitivity estimates and uncertainty extensions must be labelled according to the validation actually performed.

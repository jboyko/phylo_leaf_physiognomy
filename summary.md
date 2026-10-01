# Rosas-Puchuri et al. (2024): relevance to fossil MAP prediction

Source: *Non-linear phylogenetic regression using regularised kernels*, Methods in Ecology and Evolution 15:1611–1623, DOI [10.1111/2041-210X.14385](https://doi.org/10.1111/2041-210X.14385). Read the user-supplied 13-page PDF. Also inspected the authors' [software README](https://github.com/ulises-rosas/phylokrr).

## Central claim and implementation

The paper extends phylogenetic regression with kernel ridge regression (phyloKRR). It first transforms the predictors and response using a symmetric inverse square root of a modelled phylogenetic covariance matrix: X* = P X and y* = P y, where P'P = V^-1. It then fits regularized kernel regression in that transformed space. An RBF kernel permits nonlinear relationships and interactions; a linear kernel provides a simpler comparator. Ridge regularization penalizes complexity, and internal cross-validation chooses hyperparameters.

This is more specific than merely including taxonomy as an extra predictor. The phylogenetic covariance changes the data on which the kernel operates. The paper's regularization lambda is not our Pagel lambda; those would need distinct names in any implementation.

## Evidence and structure

The mathematical development is followed by simulations, two empirical applications, prediction-error comparisons, and feature-importance/partial-dependence analyses. Simulations use 100 trees with 500 taxa, three Brownian traits, and responses constructed to be nonlinear in the phylogenetically transformed predictors. This directly tests the method's intended transformed-space relationships. Performance generally improves over PGLS when these relationships are nonlinear; approximately linear scenarios provide a control.

The empirical examples comprise 323 fish species (biogeography and diversification) and 90 primate species (morphology). The reported improvement of about 20% in RMSE is supported by the fish example. The primate example shows no significant difference among the tested methods. The paper uses 60/40 training/test partitions and repeated comparisons, with inner folds for hyperparameter tuning. The results support the feasibility of nonlinear phylogenetic regression; they do not establish improvements for leaf climate inference or fossil transfer.

The discussion covers covariance/tree misspecification, computational scaling, alternative kernel combinations, and interpretation. Its multi-output extension predicts responses separately, so it does not itself demonstrate borrowing information from MAT to improve MAP. It does not establish species-to-site climate uncertainty coverage or an age-aware fossil prediction workflow.

## Relevance to this project

Our recent GAMs and random forests were trait-only candidates. They did not test a nonlinear trait model coupled to the phylogenetic structure that drives PIP's advantage. Their weaker performance therefore does not settle whether nonlinear phylogenetic prediction can help MAP. The paper motivates that next comparison.

However, phyloKRR's transformed-space approach is not automatically a replacement for PIP. PIP explicitly uses cross-covariance between a new occurrence and calibration taxa to adjust its prediction. A nonlinear fit accounting for phylogenetic dependence during training need not provide that same prediction-time operation.

There is also a validation/transfer issue to resolve. The software README transforms the complete X and y with `weight_data` and then splits transformed rows. A dense P means a transformed training response can depend on original responses assigned to a test set. That example is not a suitable template for claiming held-out raw site-climate prediction in this project. This observation concerns the displayed workflow and our intended target; it is not a comprehensive audit of every experiment in the publication. A new fossil has no observed climate to enter such a joint transformation.

Nonlinearity and whitening also do not generally commute: f(PX) is not P f(X). The paper's transformed-space kernel model should be distinguished from a model with nonlinear effects of original leaf traits and phylogenetically correlated residuals.

## Proposed adaptation, not an implementation of the published algorithm

A suitable project-specific model would fit a regularized nonlinear function of original, training-scaled fossil-measurable traits together with phylogenetic residual covariance. One formulation is penalized GLS over a kernel function, minimizing (y - f(X))' V^-1 (y - f(X)) plus a kernel complexity penalty, with an explicit intercept. An equivalent Gaussian-process formulation can combine a trait kernel, phylogenetic covariance, and independent residual variance. This is a proposed modelling choice, not a claim that it is exactly the paper's implementation.

For a new fossil, evaluate the trait function and use a coherently derived prediction-time phylogenetic adjustment based on its dated placement. All covariance blocks must retain the same scaffold root. Species-by-site predictions would still be averaged on the log(MAP) scale, preserving the existing fossil-local aggregation rule. Fit and tune every transformation, kernel parameter, and covariance parameter within training sites only.

The initial benchmark should compare a linear kernel and an RBF kernel under identical nested site-grouped folds, alongside current PIP and nested recalibrated PIP. The target to beat is the recalibrated benchmark (RMSE 0.422 log units), not only uncorrected PIP (0.529). Check raw site-scale error and prediction spread; if promising, then test coarsened taxonomic placement. Kernel flexibility could still yield compressed predictions, and it does not by itself fix the species-average versus local-site target distinction or interval coverage.

No kernel model has been fitted in this reading step.

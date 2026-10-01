# Nonlinear phylogenetic MAP experiment

This is an exploratory point-prediction experiment under the [prespecified
protocol](nonlinear_phylogenetic_prediction_protocol.md). It leaves the
production PIP models, fossil predictions, and uncertainty tables unchanged.
The kernel method here is an adaptation motivated by Rosas-Puchuri et al.
(2024), not their phyloKRR algorithm.

## Model and validation

The response is the natural log of species-aggregated MAP in centimetres. The
candidate model has an unpenalized intercept or linear trait mean, a linear or
RBF trait kernel, and Pagel-type residual phylogenetic covariance. All
biological covariance blocks come from the rooted extant scaffold and are
divided by one fixed scaffold root height. Numerical Cholesky jitter is kept
separate. The prediction cross block combines trait-kernel similarity with
phylogenetic covariance; no second PIP correction is added.

The existing 92 sites are evaluated in the same ten outer folds as
`code/03f_map_model_experiment.R`. Five recorded inner site folds per outer
training set select hyperparameters by pooled held-site squared error. Every
inner fit reconstructs species means, fossil-trait eligibility, trait-only
bagged imputation, and trait scaling from its own training sites. Query traits
are aggregated only within species and held site. Occurrence predictions are
averaged on the log scale before exponentiation. The initial fixed grid is
exactly the one in the protocol. Boundary selections are recorded.

The controls are the archived production PIP and its nested recalibrated
version, plus a matched-preprocessing linear-mean phylogenetic model with
`eta=0`. The production PIP's training imputer includes log(MAP), whereas the
new models use a trait-only imputer. In the first outer fold, this changes 403
of 17,352 completed training trait cells, with maximum absolute difference
0.00395. The matched control helps keep that preprocessing effect distinct
from kernel effects.

## Reproduce

From the repository root, with the pinned local `dilp` fork available:

```sh
Rscript tests/test_map_kernel_model.R
DILP_SOURCE=/path/to/dilp Rscript tests/test_map_kernel_pip_parity.R
DILP_SOURCE=/path/to/dilp PIP_KERNEL_WORKERS=3 Rscript code/03h_map_kernel_cv.R
Rscript code/03i_map_kernel_results.R
DILP_SOURCE=/path/to/dilp Rscript code/03k_map_kernel_fit_audit.R
```

`PIP_KERNEL_FOLDS=1,2` runs selected outer folds; checkpoints in
`models/map_kernel_experiment/outer_XX.rds` make these batches composable.
The result script requires all ten checkpoints and audits membership,
selection, aggregation, and comparator parity before writing tables and a
common-axis observed-versus-predicted plot. The fit audit reports Cholesky
jitter, mean columns removed for rank, and predictors removed for constancy.

After interpreting the extant comparison, the optional transfer and fossil
sensitivity scripts are:

```sh
DILP_SOURCE=/path/to/dilp Rscript code/03j_map_kernel_coarsened_transfer.R
DILP_SOURCE=/path/to/dilp PIP_KERNEL_FOSSIL_CANDIDATE=linear_mean_rbf Rscript code/04d_map_kernel_fossil_sensitivity.R
```

The fossil script requires an explicit candidate name. It pools the recorded
inner site losses to tune that family on all extant data, fits all extant
species, and predicts the two taxonomy scenarios with each fossil occurrence's
local traits and age. These are sensitivity estimates, not validated fossil
intervals. The coarsened transfer test uses only training species as taxonomic
placement anchors and evaluates genus, family, and order resolution within
the same held-site folds. It places each extant query independently at
0.00001 Ma because the scaffold's extant tip depths differ by a few numerical
branch units, causing exact age-zero grafting to fail. This 10-year offset is
negligible on the tree's million-year scale. The test does not measure
deep-time placement error.

## Results

All 92 sites have one paired held-out prediction per candidate. The numerical
and membership audits pass. RMSE and MAE below are in natural-log centimetres.
See the [common-axis observed-versus-predicted plot](../plots/map_kernel_experiment_comparison.png)
and `tables/map_kernel_experiment/site_predictions.csv` for every paired site.

| Candidate | RMSE | MAE | Mean signed error | Correlation | Predicted / observed SD |
| --- | ---: | ---: | ---: | ---: | ---: |
| Nested recalibrated PIP | 0.4220 | 0.3151 | +0.0161 | 0.753 | 0.707 |
| Linear mean + RBF + phylogeny | 0.4408 | 0.3244 | +0.1003 | 0.787 | 0.522 |
| Linear mean + phylogeny, matched preprocessing | 0.4467 | 0.3350 | +0.0966 | 0.806 | 0.468 |
| Linear kernel + phylogeny | 0.4467 | 0.3350 | +0.0971 | 0.806 | 0.467 |
| RBF kernel + phylogeny | 0.4478 | 0.3288 | +0.1068 | 0.799 | 0.482 |
| Existing PIP | 0.5295 | 0.4129 | +0.1576 | 0.750 | 0.318 |

The linear-mean RBF candidate reduces RMSE by 0.0059 (1.3%) relative to its
matched `eta=0` control. The intercept-only RBF is slightly worse than the
matched linear-kernel candidate. Neither nonlinear candidate beats nested
recalibrated PIP; the best kernel is 0.0188 (4.4%) higher in RMSE. The new
models improve substantially over the original PIP, but the matched control
shows that most of this gain is not isolated evidence for nonlinearity.
Moreover, predicted MAP remains compressed relative to observations.

Every candidate selected `lambda_phy=0.99` in all ten outer folds, the upper
edge of the prespecified grid. The linear-mean RBF selected the widest tested
bandwidth (`ell_factor=2`) ten times and `eta=10` nine times. These boundary
selections limit claims about an optimum. A wider search would be a separate
exploratory run, with its grid fixed before examining its outer errors.
No numerical jitter was needed for any of the 40 selected outer fits, and
refitting them reproduces saved occurrence predictions. All 92 site predictions
equal their saved occurrence means on the log scale. The production PIP limit
test matches existing site predictions within 1e-6 with its original
response-inclusive imputation. The real training/query joint covariance check
is positive semidefinite before numerical jitter.

The average prediction components of the linear-mean RBF are 7.56 from the
unpenalized mean, -2.43 from the trait kernel, and -0.12 from phylogeny. The
large opposing mean and kernel terms make their separate biological
interpretation unstable; the sum is the relevant point prediction.

The best nonlinear family was applied to fossils as a sensitivity estimate,
with all-extant hyperparameters chosen from pooled training-side inner losses:
`lambda_phy=0.99`, `eta=10`, and RBF bandwidth twice the extant median trait
distance. Formal-taxonomy MAP estimates range from 126.0 to 204.8 cm across
the ten fossil sites. Eight of 360 formal-only occurrences have at least one
scaled trait outside the extant training range. Two have maximum RBF
similarity below 0.5 to every calibration species; this is a diagnostic
threshold, not a validated exclusion rule. The sensitivity estimates are
in `tables/map_kernel_experiment/fossil_site_sensitivity.csv` and
`fossil_occurrence_sensitivity.csv`; the same tables include production PIP and
recalibrated PIP comparators. This does not establish which model transfers
better to fossils or provide calibrated intervals.

The coarsened held-site transfer test predicts the same 92 extant sites with
placement evidence limited to the outer training species. For the
linear-mean RBF candidate, RMSE rises from 0.4408 with the existing exact
species relationships to 0.4741 at requested genus resolution, 0.5154 at
family, and 0.5719 at order. Its matched `eta=0` control scores 0.4997,
0.5458, and 0.6122 at those resolutions. The full comparison and placement
log are in `coarsened_comparison.csv` and `coarsened_placement_log.csv`.
Of 2,682 held occurrences, 2,193 could be placed at genus level in the
genus-requested scenario; 424 fell back to family, 33 to order, and 32 to
root. This test shows material sensitivity to taxonomic resolution even on
extant taxa. It cannot validate the deep-time fossil placements.

The model-family comparison is exploratory because the winning family was
identified from these same outer scores. These folds test extant site transfer
with known relationships, while actual fossils have coarse taxonomy, extinct
lineages, and deep-time ages. The coarsened-placement check addresses only
part of that gap.

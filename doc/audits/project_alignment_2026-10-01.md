# Project and manuscript alignment audit

1 October 2026. Repository snapshot: `9c583fb`; working tree was clean at the start.

The analysis is ahead of the manuscript. The core result survives: the production PIP model improves held-out extant site predictions over the current practical site-regression and refitted DiLP baselines. The manuscript mixes older implementations and results, while the documentation contains several different versions of the project. Updating the manuscript from the Markdown files alone would carry some errors forward.

This is an audit, as requested before rewriting. The supplied ODT, existing Markdown and Word files, code, model objects, and prediction tables were not edited. Document comments were treated as historical context, not instructions. The existing `summary.md` concerns the kernel-regression paper and was preserved. No analysis was refitted and no model was promoted to production.

## What is currently implemented

| Item | Current state | Evidence |
| --- | --- | --- |
| Primary scientific application | Fossil MAT and MAP reconstruction; extant validation supplies the benchmark | [Climate prediction script](../../code/04_fossil_predictions.R) |
| Predictors | Twelve eligible fossil-measurable traits; active predictors are filtered using training data | [Authoritative trait list](../../code/01_nophy_regression.R), [validation implementation](../../code/03_loso_cv.R) |
| Extant climate sites | 92 after combining Yasuni ridgetop and upper slope | [Site grouping](../../code/site_grouping.R), [saved predictions](../../tables/loso_cv_site_predictions.csv) |
| Extant fitting unit | 1,740 operational taxon labels in the full climate calibration; grand means across training sites | [Calibration data](../../data/data_species.csv) |
| Main climate validation | Ten folds, grouping whole sites; MAT ranks assigned round-robin | [Validation implementation](../../code/03_loso_cv.R) |
| Separate LOSO experiment | 92 fits, one excluded site per fit; conditional-interval coverage check | [LOSO report](leave_one_site_out_coverage.md) |
| Fossil climate input | June 2026 source; 360 species-by-site occurrences, 346 species labels, ten analytical sites | [Fossil input](../../data/fossil_traits.csv) |
| Fossil occurrence handling | Local traits, local site age, unique occurrence tip; fixed original extant anchors | [Placement helper](../../code/fossil_placement.R), [prediction script](../../code/04_fossil_predictions.R) |
| Fossil taxonomy | Formal ranks only in the primary analysis; quoted ranks provisionally included in a separate sensitivity case | [Taxonomy rules](../../code/fossil_taxonomy.R) |
| Missing fossil climate predictors | Imputed from extant calibration; production uses the imputed PIP variant for both targets | [Fitted components](../../code/02_phy_regression.R), [prediction script](../../code/04_fossil_predictions.R) |
| Climate covariance | Calibration, CV, and fossil blocks retain the same climate scaffold root | [Root invariant test](../../tests/test_climate_covariance_root.R) |
| Site estimators | Arithmetic mean of occurrence MAT; mean occurrence ln(MAP), then exponentiate for geometric MAP in cm | [Aggregation helper](../../code/site_prediction.R), [fossil prediction](../../code/04_fossil_predictions.R) |
| LMA fossil application | Single log10(PW²/A) predictor; 192 eligible occurrences at nine sites; occurrence-specific placement; arithmetic site mean after back-transformation | [LMA prediction](../../code/04b_lma_fossil_predictions.R), [species outputs](../../tables/lma_fossil_species_predictions.csv) |
| Later MAP work | Nested recalibration, GAM/RF comparisons, kernel models, coarsened placement, and fossil sensitivity outputs exist separately | [MAP experiment](../map_model_experiment.md), [kernel experiment](../map_kernel_experiment.md) |

Species can occur at both training and held-out extant sites. Main validation measures transfer to new sites with known extant relationships. It does not directly validate extinct lineages, fossil ages, or uncertain fossil placement.

## What the current results support

The production imputed PIP model has MAT RMSE **3.4515 °C** and ln(MAP) RMSE **0.5295**, each across all 92 sites. Within the original PIP variants, imputation performs best for MAT; complete-case training performs best for ln(MAP), with RMSE **0.5122**. Both PIP variants impute missing held-out predictors and evaluate all 92 sites. Using imputation for both fossil targets is the retained production choice, not the lowest-RMSE variant for both targets.

The direct three-model comparison uses matching sites within each target:

| Model | MAT RMSE, 92 common sites | ln(MAP) RMSE, 90 common sites |
| --- | ---: | ---: |
| PIP, imputed | 3.4515 | 0.5353 |
| Refitted DiLP regression | 3.8912 | 0.5755 |
| Twelve-trait site regression, zero-filled and imputed | 4.2694 | 0.6952 |

Source: [matching-site comparison](../../tables/dana_cv_comparison.csv). MAP DiLP cannot evaluate two sites with missing required predictors. Its 90-site comparison must not be presented beside the all-92-site PIP score as though the denominators matched. The untoothed-excluded site LM is another, stronger non-phylogenetic comparator: 3.7669 °C and 0.5936 ln(MAP), across 92 sites. It is distinct from the published DiLP equations.

PGLS and PIP share the fitted regression. PGLS uses the trait component alone; PIP adds the covariance-weighted residual correction. Their comparison supports the value of that prediction-time correction. Comparing PIP with a site-trained LM changes fitting level and aggregation as well as phylogenetic treatment.

The later 92-site MAP experiments give:

| Candidate | ln(MAP) RMSE | Status |
| --- | ---: | --- |
| Nested recalibrated PIP | 0.4220 | Best observed exploratory candidate |
| Linear mean + RBF kernel + phylogeny | 0.4408 | Exploratory; small gain over matched linear control |
| Matched linear-mean phylogenetic control | 0.4467 | Exploratory control |
| Production imputed PIP | 0.5295 | Retained production estimator |

Sources: [recalibration scores](../../tables/map_model_experiment/model_comparison.csv), [kernel scores](../../tables/map_kernel_experiment/model_comparison.csv). Hyperparameters/recalibration were selected within outer training sets, but selecting a winning model family after inspecting the outer scores remains exploratory. The best tested nonlinear kernel does not beat recalibrated PIP. Its RMSE rises from 0.4408 with exact extant relationships to 0.4741, 0.5154, and 0.5719 under requested genus, family, and order placement. That coarsening test addresses taxonomic resolution, not deep-time fossil transfer.

The current univariate LMA models compare on 107 evaluable sites: PIP RMSE **0.1365**, LM **0.1468**, and PGLS **0.1473**, in log10 units. The code also evaluates multivariate LMA variants. Their smaller complete-case sample prevents direct comparison of unmatched headline scores; [the same-site table](../../tables/lma_cv_rmse_same_site.csv) supplies a 62-site comparison. The draft's description of only three tested LMA models is incomplete unless the manuscript explicitly limits reporting to the retained univariate fossil estimator.

## Manuscript map and drift

Audited source: [Phylogenetic Digital Leaf Physiognomy Manuscript.odt](</Users/jboyko/Downloads/Phylogenetic Digital Leaf Physiognomy Manuscript.odt>). SHA256: `f1f9e7204c54451ca31d4fc3845f8deb9e682a3b88b42c5cfeec79173f611356`.

The introduction motivates combining leaf physiognomy with taxonomic information and introduces climate and LMA applications. The methods explain PIP, alternative models, site validation, and fossil placement. Results compare extant errors and apply the estimator to fossils. The discussion interprets residual borrowing and begins the fossil interpretation. S1 provides the numerical example. The abstract is a placeholder, the final introduction paragraph is unfinished, fossil discussion is largely a note, and the conclusion is empty. The document therefore needs substantive completion after factual reconciliation.

| Location | Draft claim or omission | Current state / required correction |
| --- | --- | --- |
| Leaf data / fossil taxonomy | Quoted names are “not used in phylogenetic analyses” | Censor informal ranks as placement evidence; retain the occurrence using other formal ranks or root fallback. All 360 climate occurrences are retained. |
| PIP explanation | `V_inv e` is described as GLS whitening and decorrelation | It is inverse-covariance weighting. Multiplying by V⁻¹ does not by itself produce independent unit-variance residuals. The Markdown explanation already corrects this. |
| PIP explanation / empirical methods | Back-transform MAP predictions before site averaging | Average ln(MAP) predictions first, then exponentiate. This changes the site estimator and must agree with validation. |
| PIP explanation / empirical methods | The same fossil species shares one phylogenetic adjustment across sites | Every species-by-site occurrence is placed at its own age. Cross-covariance and adjustment can differ across occurrences. |
| PIP explanation | One species row is presented as required by PGLS | Grand-mean calibration is a choice of this implementation. It loses within-taxon site variation; unnamed labels are not all verified species. |
| Validation design | Snake-pattern fold assignment | Actual code assigns MAT ranks modulo ten, round-robin. The `loso` filename is historical; default validation is ten-fold. |
| Validation / missing data | Imputation is learned on the full dataset before CV | All validation imputers are learned inside their training folds. Production PGLS training imputation includes its known training response; held-out imputation uses training-fitted trait-only models. Avoid claiming all imputers are identical. |
| Validation / complete cases | Complete-case filtering uniformly removes held-out sites | PGLS/PIP CC means complete training rows, with imputed held-out traits. They retain 92 sites. LM CC variants omit incomplete held-out inputs and have smaller evaluated sets. |
| Validation / model configurations | Specimen site LM directly weights all specimens and differs from sp+zero | Both currently average processed morphotype means and are duplicate configurations. They provide one comparison, not two independent methods. |
| Validation / covariance | Full pruned extant tree supplies the covariance | Climate covariance is subset from the rooted scaffold matrix. Recomputing VCV after pruning discards shared stem history and is incompatible with current fossil blocks. |
| Climate Results / Tables 1–3 | Older errors, slopes, and 93-site captions | Replace the numerical tables together from current outputs. Some draft tables also omit model rows despite captions referring to all twelve configurations. Include the refitted DiLP baseline deliberately. |
| LMA Results / Tables 4–5 | PIP 0.130, LM 0.146, PGLS 0.144; all 108 sites | Current retained univariate comparison is 0.1365 / 0.1468 / 0.1473 on 107 sites. Distinguish the foldable data from each estimator's evaluable sites. |
| Empirical methods / Tables 6–7 | 361 climate rows; eleven fossil sites; PL1 and PL2 separate; older ages | Current climate input has 360 occurrences at ten sites. Combine Palacio before aggregation at 64.08 Ma. Replace ages, predictions, and counts together. LMA has 192 occurrences at nine sites; Cerrejon lacks the usable fossil petiole predictor. |
| Empirical Results | Old fossil climate/LMA predictions | Current formal-only climate MAT ranges 14.1–21.6 °C and MAP 146–216 cm in the rounded table. Current LMA ranges 74.7–125.9 g m⁻². These remain predictions, not demonstrations of fossil accuracy. |
| Empirical uncertainty | CV RMSE is treated as an expected fossil error, with species count implying precision | RMSE is an extant performance benchmark. Fossil uncertainty depends on joint errors, ages, placement, and input uncertainty; species count alone is insufficient. Conditional site intervals fail coverage. |
| Reproduction statement | Repository and `dilp_pgls()` are interchangeable reproductions | The available local package has a different API from the parity test, and the pinned public commit does not contain the revised repository prediction implementation. Establish parity/versioning before making this claim. |
| Discussion / compression | A small correction is said to pull a prediction toward the training mean | A small correction leaves the PGLS trait prediction largely unchanged. Compression is observed, but this proposed causal explanation is not established by the algebra or error decomposition. |
| Discussion / mechanism | Similar LMA bias is said to establish recovery of biological signal beyond regularization | Better held-out error supports predictive utility. Error decomposition does not identify a biological mechanism or rule out regularization. Leaf habit is a hypothesis, not a demonstrated recovered variable. |
| Discussion / sparse collections | Correction is claimed to provide benefits with few species or uncertain taxonomy | Specimen count does not index the correction, but that fact alone does not establish accuracy or interval coverage for sparse fossil assemblages. |
| Main text overall | Later uncertainty, recalibration, and kernel work absent | Decide which work belongs in the paper, supplement, or research record before incorporating it. A production-method description should not silently inherit experimental estimators. |

### The worked example

Keep its explanatory structure. It separates a 17.5 °C trait prediction from a negative residual correction and arrives at about 15.3 °C. The assumed regression coefficients are explicitly illustrative; with two observations and two freely fitted coefficients, a newly fitted regression would interpolate the observations and have no residual correction.

There is a small rounding issue: using the exact inverse of the stated covariance gives a correction of −2.2143 °C and a final value of 15.2857 °C. Using the displayed rounded inverse gives −2.203 °C and 15.297 °C. Both round to 15.3 °C. Resolve the intermediate precision during editing without changing the example's point. Use “covariance-weighted residuals” consistently; both extant residuals and their joint covariance contribute to the correction.

## The uncertainty problem is still open

The joint conditional-error calculation is implemented and passes generated-output checks, including full covariance propagation within sites and coefficient uncertainty. Numerical correctness does not establish coverage for actual site climate.

| Nominal 95% interval | MAT covered / 92 | ln(MAP) covered / 92 |
| --- | ---: | ---: |
| Conditional, ten-fold | 35 (38.0%) | 29 (31.5%) |
| Conditional, leave one site out | 35 (38.0%) | 25 (27.2%) |
| Shared-discrepancy augmentation, reported fold-exclusion check | 89 (96.7%) | 86 (93.5%) |

Sources: [LOSO comparison](../../tables/leave_one_site_out/coverage_comparison.csv), [augmented-interval summary](../../tables/pip_cv_calibrated_coverage.csv). The audit independently recounted the ten-fold conditional coverage from numeric endpoints.

The third row needs a qualification missing from the current prose. [The diagnostic script](../../code/03d_uncertainty_diagnostics.R) estimates the extra variance from other folds' saved out-of-fold errors. It excludes the evaluation fold's residuals directly, but the models producing those other residuals were trained using the evaluation fold's sites. Consequently, the full interval-calibration procedure is not nested entirely within the outer training set. This is a methodological inference from the fold construction; the reported counts themselves are correct. Its effect on coverage has not been measured here.

Before claiming independent coverage of the augmented estimator, generate its calibration residuals through inner validation restricted to each outer training set, then evaluate untouched outer sites. The later nested MAP recalibration experiment does follow that training restriction, but it evaluates point predictions and does not validate these intervals.

The uncertainty document also moves from “do not replace fossil error bars” to “report the calibrated interval.” That is an unresolved recommendation, not a settled result. Calling the augmented fossil intervals a “lower bound on true uncertainty” is not an established bound. They omit uncertainties and lack fossil-transfer validation, but omission alone does not prove their widths bound the true error distribution. Describe them as experimental until the validation target and procedure are agreed.

## Where the documentation has drifted

| File or group | Current problem | Treatment in a later alignment pass |
| --- | --- | --- |
| [AGENTS.md](../../AGENTS.md) | Hardcoded-path instructions, routine LOSO description, legacy package loading, and imputation recommended as best for both targets | Preserve collaborator constraints; refresh operational descriptions from current code and results. |
| [CLAUDE.md](../../CLAUDE.md) | More recent CV description, but old 93-site scores and fossil comparators that were removed | Give it the same factual contract as AGENTS.md. |
| [README](../../README.md) | Main climate explanation mostly current; LMA package-loading statement and site-aggregation axis labels remain stale | Keep it an execution/data guide; link manuscript facts rather than duplicating historical scores. |
| [Domain context](../../docs/agents/domain.md) | Says fossil MAP averages back-transformed predictions | Correct to log-scale averaging when prose edits are authorized. |
| [PIP explanation](../fossil_predict.md) | Closest to current point-prediction explanation; worked example worth retaining | Update facts only where needed; do not presume this establishes the desired manuscript voice. |
| [Validation](../model_validation.md) | Current point-score tables, but pruned-tree climate covariance prose and incomplete LMA/model-count explanation | Reconcile methods with the exact retained estimators. |
| [Empirical application](../empirical.md) | Pre-root-correction point values; wrong formal placement totals; pruned-tree climate covariance; package equivalence overstated | Refresh from current outputs after decisions below. |
| [Discussion](../discussion.md) | Old RMSE and in-sample DiLP reference; residuals called decorrelated | Update evidence and distinguish hypotheses from measured outcomes. |
| [Uncertainty](../prediction_uncertainty.md) | Conflicting reporting recommendations, overstated nesting, causal certainty, and unsupported lower-bound language | Agree validation/reporting status first. |
| [Release README](../../release/README.md) | States arithmetic MAP aggregation, omits scaffold inputs from some dependency rows, and says LMA lacks occurrence updates | Reconcile with release code rather than relying on historical descriptions. |
| Dated Dana reports and audits | Different snapshots; some contain later edits mixed into older reports | Preserve as history, label what was superseded, and link one current status source. |
| Word extracts | Separate tracked copies; September climate trial explicitly did not refresh them | Regenerate only after their editable sources and manuscript scope agree. Their contents were not fully audited here. |
| `summary.md` | Reading notes for Rosas-Puchuri et al., not a project summary | Preserve the notes; give their purpose an explicit label/link. |

For example, `doc/empirical.md` reports formal-only placements as 192 root, 99 family, 49 genus, and 20 order. The current saved log gives **195 root, 105 family, 37 genus, and 23 order**, totaling 360. Root attachment is successful prediction fallback with weak taxonomic resolution, not precise placement. The older counts must not be reused.

All 37 shared main/release R implementations match as parsed expressions after accounting for three known LMA bootstrap-loader differences. The four bytewise differences consist of those loaders, comments, and whitespace; they do not establish a separate modelling implementation. The release documentation and package status still prevent treating all entry points as interchangeable.

## What needs a decision, rather than a copy edit

1. **Paper scope and estimator.** Retain production PIP as the main method unless explicitly deciding to adopt a later candidate. Decide whether MAP recalibration, kernels, and coarsened transfer belong in the main paper or supplement. Their existing fossil sensitivity tables do not establish fossil accuracy.
2. **Uncertainty reporting.** Conditional intervals fail the site-climate target. Augmented intervals need properly nested calibration validation. Choose which experimental diagnostics to report and what uncertainty claim the paper can defend.
3. **Unnamed calibration labels.** The review inventory contains 213 labels without species epithets; 79 pool more than one site. Retaining these pools is the current convention. Changing biological identity assumptions requires collaborator review and a full rebuild/refit, not wording changes.
4. **LMA presentation and estimand.** Climate uses a rooted scaffold; LMA uses its own pruned calibration tree. The saved univariate LMA training covariance agrees with that LMA prediction-base tree to about 1e-13, so this audit does not establish the former climate root mismatch in LMA. However, LMA CV scores mean log10 predictions, while fossil site outputs average after back-transformation. Explain this distinction before attaching log-scale CV errors to arithmetic site means. Decide whether to report only the retained univariate fossil method or also the multivariate experiments.
5. **Reproduction route.** Identify the repository/release commit that reproduces the selected estimator. Bring `dilp_pgls()`, exported model objects, and tests into parity before claiming package reproduction. The local package currently accepts only `specimen_data`; the repository parity test passes `taxonomy_scenario`.

These decisions can be made without reopening settled fossil trait eligibility, site grouping, units, occurrence ages, or formal-only primary taxonomy.

## Verification performed for this audit

Passed: climate covariance-root test; real fossil placement and insertion-order checks; site-prediction aggregation; fossil taxonomy; 92-site grouping/revised ages; generated CV/fossil uncertainty outputs. Independent Python checks reproduced the conditional coverage counts and all ten unrounded fossil MAT/log-MAP site means against the saved uncertainty estimates.

The climate occurrence integration test **failed before prediction**, because its isolated temporary project does not provide the now-required `tables/pip_cv_site_uncertainty.csv`. The production outputs exist; the test fixture has drifted after addition of site-discrepancy loading. This failure is not evidence of wrong fossil point estimates, but it prevents claiming that the current integration test passes unchanged.

The LMA occurrence integration test passed in its isolated temporary project. No full pipeline rerun, package publication, external literature verification, or exhaustive export/layout audit was performed.

## Order of the alignment work

First agree the estimator, uncertainty claim, and treatment of later experiments. Then freeze one factual record with source tables, denominators, site estimators, and model versions. Reconcile guidance and execution docs against that record. Refresh manuscript tables and dependent claims together. Finally edit the main text paragraph by paragraph, retaining the worked example and choosing the writing voice with James rather than inheriting it from old Markdown. Regenerate Word extracts and the release documentation after the substantive text agrees.

Most factual drift can be corrected without rerunning the production climate analysis. Changing calibration identities, adopting a new estimator, or establishing augmented-interval coverage requires analytical work. Those should not be concealed inside a manuscript rewrite.

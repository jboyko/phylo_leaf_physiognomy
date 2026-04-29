# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview and Primary Goal

This project predicts paleoclimate (MAT = mean annual temperature, MAP = mean annual precipitation) from fossil leaf morphology using phylogenetically-informed prediction (PIP). **The primary deliverable is fossil climate prediction** — not the extant species model comparison. Every modelling decision should be evaluated against the question: "can this be applied to a fossil leaf?"

The extant species pipeline (scripts 00–03) is a calibration and validation benchmark. Script 04 is the actual scientific product.

---

## Key Decisions from Correspondence with Dana Royer

Dana Royer is the fossil leaf expert and domain collaborator. The following decisions were made explicitly in correspondence with him and must not be silently reversed:

### 1. Fossil-measurable traits only
**All models are restricted to traits that can be measured from fossil leaf specimens.** Dana provided the definitive list (pers. comm.). This means `evergreen` — despite being the single strongest MAT predictor in early runs — is excluded. It cannot be determined from a fossil. The same applies to petiole-based traits (`petiole.width.cm`, `petiole.area.cm2`), `major.length.cm`, `compact`, `shape.factor`, `perim.area.cm2`, and raw tooth counts.

Training on all available extant traits and then imputing unmeasurable ones for fossils would silently undermine the fossil application. The trait set must be the same for training and prediction.

**The 12 fossil-measurable predictors (defined in `01_nophy_regression.R`):**
- Always measurable: `pw2.a.ratio`, `ln.leaf.area.mm2`, `feret.diam.ratio`, `margin.score`
- Toothed leaves (0 or 1 for untoothed): `perim.ratio`, `teeth.perimeter.percm`, `teeth.interior.percm`, `avt.tooth.area`, `tooth.area.blade.area.ratio`, `tooth.area.perimeter`, `tooth.area.interior`, `teeth.blade.area.ratio`

Note: `ln.leaf.area.mm2` represents Dana's `leaf.area.cm2` — same measurement, log-scaled for regression.

### 2. Tooth trait filling
For confirmed untoothed leaves (cell blank **AND** `margin.score == 1`), tooth trait NAs are biologically real zeros — set them explicitly **before** aggregating to species means in `00_data_cleaning.R`. Without this, tooth traits appear as ~67% NA and get excluded from all models. Dana confirmed:
- Set to **0**: `primary.teeth.number`, `secondary.teeth.number`, `teeth.number`, `teeth.perimeter.percm`, `teeth.interior.percm`, `tooth.area.cm2`, `avt.tooth.area`, `tooth.area.perimeter`, `tooth.area.interior`, `tooth.area.blade.area.ratio`, `teeth.blade.area.ratio`
- Set to **1**: `perim.ratio` (untoothed leaves have a smooth perimeter, ratio = 1)

Both conditions must be met: the cell is blank AND `margin.score == 1`. This is a security measure Dana requested for when we apply the model to fossil data.

### 3. Species-level analysis
We use species-level trait means (not site means) as the basis for modelling, because we are adding phylogeny and phylogeny operates at the species level. Site-level predictions can be derived afterwards by aggregating species-level predictions within a site. Dana agreed this is the right approach.

### 4. Climate targets
MAT and log(MAP) are the primary targets. Other climate variables (coldest month temperature, growing degree days, etc.) covary with MAT/MAP and are not modelled separately.

### 5. Fossil placement
Fossils are commonly known only to family or order level, or may belong to extinct genera. The code handles this via a genus → family → order MRCA fallback, with root placement as the final fallback. Dana confirmed this is appropriate. Fossil age is used to set edge length so the fossil tip sits at the correct time depth.

When a fossil predates the crown age of its placement clade, the `PLACEMENT_FALLBACK` flag in `04_fossil_predictions.R` controls behaviour: `"ancestral_branch"` (default) walks up the tree to find the branch alive at the fossil's age and splits it there; `"node"` attaches at the MRCA with a minimal edge (0.001 Ma). The ancestral branch approach is preferred as it respects the fossil's age estimate.

### 6. `teeth.interior.percm` in MAT model
Dana's original temperature model included "number of teeth / internal perimeter" (`teeth.interior.percm`). This was absent from early runs because tooth traits were excluded. After the tooth-filling fix and restriction to fossil-measurable traits, ElasticNet now selects it for MAT. This is consistent with Dana's prior work.

### 7. `internal.perimeter.cm`
This is a raw tooth-linked measurement (~67% NA) not on Dana's fossil list. It remains excluded. The ratio forms (`teeth.perimeter.percm`, `teeth.interior.percm`) are included instead and are more appropriate for cross-specimen comparison.

### 8. MAP units
The training data `map` column (from `RoyerLeafShapeClimateDataFixedNames_June2012.csv`) is in **centimetres**, not millimetres. The range is approximately 19–680 cm (mean ~215 cm). All model MAP outputs are therefore in cm. The published DiLP MAP estimates in Peppe et al. (2011) are also in cm and are directly comparable to model output. Do not convert.

### 9. Fossil site ages
Ages used for phylogenetic placement of Peppe et al. (2011) fossil sites (midpoints of published ranges):

| Site | Age (Ma) |
|------|----------|
| Fox Hills | 66.5 |
| Williston Basin I | 64.75 |
| Williston Basin II | 63.5 |
| Williston Basin III | 59.75 |
| Palacio de los Loros PL1 & PL2 | 61.7 |
| Cerrejon | 58.0 |
| Hubble Bubble | 55.8 |
| Laguna del Hunco | 51.9 |
| Republic | 49.4 |
| Bonanza | 47.3 |

### 10. Species-level vs site-level aggregation — a key methodological finding
`04_fossil_predictions.R` compares four prediction approaches for each fossil site:

1. **LM sp** — LM trained on extant species means; fossil species grand means predicted per species, averaged to site
2. **LM site** — LM trained on extant site means; fossil site-mean traits predicted directly per site
3. **PIP sp** — PIP using fossil species grand means, averaged to site
4. **PIP+site** — PIP using site-specific fossil species means, averaged to site (PIP sp and PIP+site share the same phylogenetic adjustment; they differ only in which trait values go into X_fossil)

**Updated finding (10-fold LOSO CV, stratified by MAT, seed = 42):** PIP site is the best
performer for both MAT and MAP. PGLS alone (without the PIP covariance correction) is
substantially worse than LM site — the correction term is what drives the improvement.

| Model                  | MAT RMSE | log(MAP) RMSE | n_sites |
|------------------------|----------|---------------|---------|
| pip_site_impute        | 3.638    | 0.542         | 92      |
| pip_site_cc            | 3.653    | 0.545         | 92      |
| lm_site_peppe_impute   | 3.913    | 0.572         | 92      |
| lm_site_peppe_cc       | 3.943    | 0.581         | 90      |
| lm_site_specimen       | 4.109    | 0.603         | 92      |
| lm_sp                  | 5.631–5.662 | 0.651–0.652 | 92   |
| pgls_site              | 6.319–6.362 | 0.674–0.675 | 92   |

PIP is acting as a regularization procedure — it encourages simpler, less overfit models and
borrows signal from phylogenetically close extant relatives. The earlier finding that "LM site
outperforms PIP" was based on a comparison against published DiLP estimates at specific fossil
sites, not a CV benchmark. The LOSO CV is the correct evaluation.

**Practical implication**: Use PIP (impute variant) as the recommended model for fossil
climate reconstruction.

---

## Running the Analysis

Scripts must be run in order. Each sets `setwd()` to a hardcoded path — update line 1 when working on a new machine.

```r
Rscript code/00_data_cleaning.R      # Fill tooth traits, aggregate, build scaffold tree
Rscript code/01_nophy_regression.R   # Fit LM/ENet/RF species-level; LM site-level (6 configs)
Rscript code/02_phy_regression.R     # Fit PGLS impute+cc variants, save pip_components.rds
Rscript code/03_loso_cv.R            # LOSO CV across all 12 model configs, RMSE tables, model summaries
Rscript code/04_fossil_predictions.R # Predict MAT/MAP for fossils via PIP
Rscript code/05_sample_size_analysis.R  # Bootstrap RMSE vs N specimens (site LM vs PIP)
```

---

## Pipeline Architecture

**`00_data_cleaning.R`** — Fills tooth trait NAs before aggregation (see section 2 above). Aggregates to species-level means and three site-level variants (see below). Builds a family-level angiosperm backbone (2 crown tips per family, ~943 tips across 515 families) from the full WCVP tree, then grafts training species onto that small scaffold. All getMRCA/bind.tip calls operate on the small backbone (~1000–3000 tips) rather than the full 123k-tip tree. `h = max(nodeHeights())` is computed once before the loop — constant because all tips reach the crown age. A `BUILD_PHYLOGENY` flag at the top (default `FALSE`) controls whether phylogeny sections run — set to `TRUE` only when the tree or species list changes. Outputs:
- `data/data_species.csv` — species-level means
- `data/dat_site.csv` — site-level means: specimen → site directly (zero-filled tooth traits; original approach)
- `data/dat_site_sp_zero.csv` — site-level means: specimen → species-within-site → site (zero-filled; species-weighted)
- `data/dat_site_peppe.csv` — site-level means: specimen → species-within-site → site, **no zero-fill** so untoothed species are excluded from tooth trait averages via `na.rm = TRUE`; matches Peppe et al. (2011)
- `data/tre_pruned.tre` — training species only; used for VCV in PGLS/PIP
- `data/tre_scaffold.tre` — family backbone + training species; used for fossil grafting in `04_`
- `data/name_table_full.csv` — genus/family/order for all 123k WCVP tips (read by `02_`, avoids reloading the full tree)

**`01_nophy_regression.R`** — Explicitly defines the 12 fossil-measurable predictors (`fossil_traits` vector at top of script). Pre-imputes once via `bagImpute` before CV for species-level models. Fits LM, ElasticNet, and RF at **species level** (saved to `models/nophy_models.rds`; ENet retained because `02_` no longer uses it for trait selection, but it remains available). Fits **LM only** at site level across all six dataset × imputation combinations (3 datasets × bagImpute/complete-case), stored in `models/site_models.rds` under `$configs`. Backward-compatible top-level keys (`$mat$LM`, `$log_map$LM`, `$impute_model`, `$pred_names`) point to the original `specimen_impute` config and are used by `05_sample_size_analysis.R`.

**`02_phy_regression.R`** — Reads `name_table_full.csv` (no need to reload the 123k tree). Uses the same `pred_names` as the LM in `01_` (loaded from `nophy_models.rds`) — **no ENet-based trait selection**, so phylogenetic and non-phylogenetic models are directly comparable. Fits PGLS for MAT and MAP via `pglmEstLambda()` in two variants: `impute` (bagImpute) and `cc` (complete-case, phylogeny pruned to complete species). Saves `models/pip_components.rds` with all variants under `$configs`; backward-compatible top-level keys point to the `impute` variant.

**Lambda transformation convention**: `V_lam` is computed to match the `lamTrans()` function used internally by `pglmEstLambda()` — off-diagonals are multiplied by λ, diagonal is left unchanged (`diag(V_lam) <- diag(phylomat)`). Do NOT use `diag(V_lam) <- diag(V_lam) + (1 - lambda)`, which adds a dimensionless nugget incompatible with the VCV scale and inconsistent with fitting.

**`03_loso_cv.R`** — Leave-one-site-out (LOSO) cross-validation across all 12 model configurations (6 model types × bagImpute/complete-case). Holds out each site entirely, retrains all models from scratch on the remaining sites, and predicts the held-out site. Stricter than analytical LOOCV — entire sites are withheld and models are retrained each fold, which better mirrors the fossil prediction setting. Column naming convention: `{method}_{train_level}_{pred_level}_{config}_{na}_{target}`. Saves per-fold model objects to `models/loso_cv_fold_XX.rds`. Outputs: `tables/loso_cv_rmse.csv`, `tables/loso_cv_site_predictions.csv`, `tables/loso_cv_model_coefs.csv`, `tables/loso_cv_model_fit.csv`.

**`04_fossil_predictions.R`** — The primary scientific output script. Loads `pip_components.rds` and `fossil_traits.csv`. Grafts fossils onto `tre_scaffold.tre` (family backbone ensures any angiosperm family can be placed). Prunes to training + fossils, computes cross-covariance V[training, fossil] via `vcv()`, applies PIP formula. Reuses `pip$V_lam_mat` and `pip$V_lam_map` — valid because covariance between two training species is independent of what other taxa are in the tree. Writes `tables/fossil_predictions.csv`.

**`05_sample_size_analysis.R`** — Bootstraps N specimens (N = 1–50) from each training site and compares RMSE for site-level LM vs site-aggregated PIP. For each replicate and site, the LM receives subsample-mean traits; PIP averages LOOCV predictions for the species represented in the subsample. Writes `tables/sample_size_results.csv` and `tables/site_specimen_counts.csv`. Key caveat: the bootstrap comparison is not apples-to-apples at large N — the LM is tested on noisier inputs than it was trained on (full site means), so it appears to never catch PIP in the bootstrap curves. The CV reference points (diamonds in the report figure) show the LM's true full-site performance, where it does beat PIP for MAT. The plot is still being refined — see pending work below.

**`code/Phylogenetically-Informed_Predictions_Source.R`** — PGLS library from Freckleton (2015), modified by Gardner et al. (2024). Key functions: `pglm()`, `pglmEstLambda()`, `pglmPredictMissing()`, `weights.p()`.

**`report.Rmd`** — R Markdown report covering the full methodology and results. Renders to `report.html` via `rmarkdown::render("report.Rmd")`. Note: requires pandoc; knit from RStudio or set `RSTUDIO_PANDOC` on the command line.

---

## Key Data Files

| File | Description |
|------|-------------|
| `data/RoyerLeafShapeClimateDataFixedNames_June2012.csv` | Raw Royer leaf morphology + climate data |
| `data/best_wcvp.tre_dated` | Dated angiosperm phylogeny (tip labels: `order_family_genus_species`) |
| `data/name_table_full.csv` | Parsed order/family/genus/species for all WCVP tips; written by `00_` |
| `data/data_species.csv` | Species-level averaged traits; output of `00_` |
| `data/tre_pruned.tre` | Training species only; used for VCV computation |
| `data/tre_scaffold.tre` | Family backbone + training species; used for fossil grafting |
| `data/fossil_traits.csv` | User-provided fossil specimen traits |
| `data/dat_site.csv` | Site-level means: specimen → site (zero-filled tooth traits); output of `00_` |
| `data/dat_site_sp_zero.csv` | Site-level means: species → site (zero-filled; species-weighted); output of `00_` |
| `data/dat_site_peppe.csv` | Site-level means: species → site, untoothed excluded from tooth averages (Peppe 2011); output of `00_` |
| `models/nophy_models.rds` | Species-level caret models (LM, ENet, RF) + `impute_model` + `pred_names`; from `01_` |
| `models/site_models.rds` | Site-level LM models across 6 configs under `$configs`; backward-compat keys + `impute_model` + `pred_names` at top level; from `01_` |
| `models/pip_components.rds` | All PIP components for fossil prediction; from `02_` |
| `tables/loso_cv_site_predictions.csv` | Per-site LOSO predictions for all 12 model configs; from `03_` |
| `tables/loso_cv_rmse.csv` | LOSO RMSE summary across all configs and targets; from `03_` |
| `tables/loso_cv_model_coefs.csv` | Per-fold coefficient estimates for all models; from `03_` |
| `tables/loso_cv_model_fit.csv` | Per-fold R², RMSE, lambda for all models; from `03_` |
| `tables/sample_size_results.csv` | Bootstrap RMSE vs N specimens for LM and PIP; from `05_` |
| `tables/site_specimen_counts.csv` | Median/mean specimens per site; from `05_` |

## R Package Dependencies

`ape`, `phytools`, `caret`, `glmnet`, `ranger`, `pdp`, `gridExtra`, `ggplot2`, `mvtnorm`, `MASS`, `rmarkdown`, `kableExtra`

## Git Commits

Do not add Claude as a co-author in commit messages.

## Pending Work

- **Sample size plot (`05_`/`report.Rmd`)**: The bootstrap comparison (filled circles) is not apples-to-apples because the LM is tested on noisier inputs than it was trained on. The CV reference diamonds show where the LM ends up at full site size, but the visual is not yet fully satisfying. A proper fix would retrain the LM within the bootstrap (leave-one-site-out CV at each N), so both methods are evaluated on held-out sites using the same N-specimen inputs. This is more expensive but would give a fair crossover curve.

## Miscellaneous

- The `old/` directory contains prior implementations and should not be modified.
- `pw2.a.ratio` is included in the fossil trait set but Dana noted "most folks don't measure it" — predictions for fossils lacking this trait rely on imputation.
- If adding new predictors, check Dana's fossil list first. The fossil-measurable constraint is a hard design requirement, not a preference.

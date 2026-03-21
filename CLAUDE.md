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

---

## Running the Analysis

Scripts must be run in order. Each sets `setwd()` to a hardcoded path — update line 1 when working on a new machine.

```r
Rscript code/00_data_cleaning.R      # Fill tooth traits, aggregate, build scaffold tree
Rscript code/01_nophy_regression.R   # Fit LM, ElasticNet, RF (fossil traits only)
Rscript code/02_phy_regression.R     # Fit PGLS (MAT + MAP), save pip_components.rds
Rscript code/03_comparison.R         # LOOCV, RMSE tables, importance, plots
Rscript code/04_fossil_predictions.R # Predict MAT/MAP for fossils via PIP
```

---

## Pipeline Architecture

**`00_data_cleaning.R`** — Fills tooth trait NAs before aggregation (see section 2 above). Builds a family-level angiosperm backbone (2 crown tips per family, ~943 tips across 515 families) from the full WCVP tree, then grafts training species onto that small scaffold. All getMRCA/bind.tip calls operate on the small backbone (~1000–3000 tips) rather than the full 123k-tip tree. `h = max(nodeHeights())` is computed once before the loop — constant because all tips reach the crown age. Outputs:
- `data/data_species.csv` — species-level means
- `data/tre_pruned.tre` — training species only; used for VCV in PGLS/PIP
- `data/tre_scaffold.tre` — family backbone + training species; used for fossil grafting in `04_`
- `data/name_table_full.csv` — genus/family/order for all 123k WCVP tips (read by `02_`, avoids reloading the full tree)

**`01_nophy_regression.R`** — Explicitly defines the 12 fossil-measurable predictors (`fossil_traits` vector at top of script). Pre-imputes once via `bagImpute` before CV. Uses `"impurity"` importance for RF (fast; `"permutation"` is slow and unnecessary for ranking). Saves `models/nophy_models.rds`. No exports or plots — those live in `03_`.

**`02_phy_regression.R`** — Reads `name_table_full.csv` (no need to reload the 123k tree). Extracts non-zero ElasticNet coefficients → active traits → PGLS formula. Fits PGLS for MAT and MAP via `pglmEstLambda()`. Saves `models/pip_components.rds` containing: beta, lambda, V_lam, residuals, design matrices, imputation models (full and traits-only), taxonomy, name_table_full, and pruned tree.

**Lambda transformation convention**: `V_lam` is computed to match the `lamTrans()` function used internally by `pglmEstLambda()` — off-diagonals are multiplied by λ, diagonal is left unchanged (`diag(V_lam) <- diag(phylomat)`). Do NOT use `diag(V_lam) <- diag(V_lam) + (1 - lambda)`, which adds a dimensionless nugget incompatible with the VCV scale and inconsistent with fitting.

**`03_comparison.R`** — All reporting in one place. Vectorised PIP LOOCV via Schur complement: `ŷ₋ᵢ = yᵢ − (Kε)ᵢ/Kᵢᵢ` where `K = V⁻¹`. One matrix solve replaces N solves — O(N³) vs naive O(N⁴). Exports: nophy CV summaries, RF variable importance, PDP plots (top 4 per target), four-model RMSE tables, scatter plots.

**`04_fossil_predictions.R`** — The primary scientific output script. Loads `pip_components.rds` and `fossil_traits.csv`. Grafts fossils onto `tre_scaffold.tre` (family backbone ensures any angiosperm family can be placed). Prunes to training + fossils, computes cross-covariance V[training, fossil] via `vcv()`, applies PIP formula. Reuses `pip$V_lam_mat` and `pip$V_lam_map` — valid because covariance between two training species is independent of what other taxa are in the tree. Writes `tables/fossil_predictions.csv`.

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
| `models/nophy_models.rds` | Caret model objects from `01_` |
| `models/pip_components.rds` | All PIP components for fossil prediction; from `02_` |

## R Package Dependencies

`ape`, `phytools`, `caret`, `glmnet`, `ranger`, `pdp`, `gridExtra`, `ggplot2`, `mvtnorm`, `MASS`, `rmarkdown`, `kableExtra`

## Git Commits

Do not add Claude as a co-author in commit messages.

## Miscellaneous

- The `old/` directory contains prior implementations and should not be modified.
- `pw2.a.ratio` is included in the fossil trait set but Dana noted "most folks don't measure it" — predictions for fossils lacking this trait rely on imputation.
- If adding new predictors, check Dana's fossil list first. The fossil-measurable constraint is a hard design requirement, not a preference.

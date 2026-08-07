# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview and Primary Goal

This project predicts paleoclimate (MAT = mean annual temperature, MAP = mean annual precipitation) from fossil leaf morphology using phylogenetically-informed prediction (PIP). **The primary deliverable is fossil climate prediction** — not the extant species model comparison. Every modelling decision should be evaluated against the question: "can this be applied to a fossil leaf?"

The extant species pipeline (scripts 00–03) is a calibration and validation benchmark. Script 04 is the actual scientific product.

## Running the Analysis

Scripts must be run in order. Each sources `code/setup.R` as its first line, which locates the repository root automatically (by walking up from the working directory looking for `README.md` alongside a `code/` directory), `setwd()`s there, and creates the gitignored `models/`, `tables/`, `plots/` output directories. There is no hardcoded path to edit — run scripts with the working directory set anywhere inside the repo, e.g. `Rscript code/00_data_cleaning.R` from the repo root.

```r
Rscript code/00_data_cleaning.R      # Fill tooth traits, aggregate, build scaffold tree
Rscript code/00c_fossil_data_cleaning.R # Build fossil traits from April 2026 leaf data
Rscript code/01_nophy_regression.R   # Fit non-phylogenetic models (species + site level)
Rscript code/02_phy_regression.R     # Fit PGLS, save pip_components.rds
Rscript code/03_loso_cv.R            # 10-fold site-grouped cross-validation (see naming note below)
Rscript code/04_fossil_predictions.R # Formal-only primary + informal sensitivity
Rscript code/05_visualizations.R     # Figures from 10-fold site-grouped CV outputs
Rscript code/06_update_dilp_package.R # Regenerate dilp package model objects + constants
```

`code/Phylogenetically-Informed_Predictions_Source.R` is the PGLS library from Freckleton (2015), modified by Gardner et al. (2024).

A parallel `*b*` / LMA pipeline (`00b`–`03b`, `04b_lma_*`, `05b_lma_*`) mirrors these steps to predict leaf mass per area instead of climate, and uses its own calibration inputs.

### `dilp` dependency

The pipeline depends on the unpublished fork `jboyko/dilp`, pinned to commit `b29be909355f7edb7809bf718f189b7a53b789d2` — the CRAN release of `dilp` carries different regression constants and will not reproduce these results. Install it with:

```r
remotes::install_github("jboyko/dilp", ref = "b29be909355f7edb7809bf718f189b7a53b789d2")
```

Every script that needs it calls `pip_require_dilp()` (defined in `code/setup.R`), which errors out with the install command above if `dilp` is not on the library path. This replaces the previous three inconsistent loading styles (`library(dilp)`, `devtools::load_all("~/dilp")`, and a `requireNamespace()`/fallback combo).

**Cross-validation naming note**: `03_loso_cv.R` (and `03b_loso_cv.R`) and their outputs are named `loso` for continuity, but the procedure is **10-fold cross-validation with sites as the grouping unit**, not leave-one-site-out. Sites are ranked by site-level MAT (or, for the LMA pipeline, site-mean log10 LMA) and dealt into 10 folds, so roughly 9 sites are held out per fold; whole sites are always held out together and every model is refit from scratch on the remaining sites. File names, columns, and log lines that reference "loso" describe this 10-fold site-grouped procedure, not true leave-one-site-out.

**Lambda transformation gotcha**: `V_lam` must match the `lamTrans()` function used internally by `pglmEstLambda()` — off-diagonals multiplied by λ, diagonal left unchanged (`diag(V_lam) <- diag(phylomat)`). Do NOT use `diag(V_lam) <- diag(V_lam) + (1 - lambda)`, which adds a dimensionless nugget incompatible with the VCV scale and inconsistent with fitting.

## Key Decisions from Correspondence with Dana Royer

Dana Royer is the fossil leaf expert and domain collaborator. These decisions were made explicitly with him and must not be silently reversed.

### Fossil-measurable traits only
**All models are restricted to traits that can be measured from fossil leaf specimens** — Dana provided the definitive list (pers. comm.). The authoritative trait set is the `fossil_traits` vector in `01_nophy_regression.R`. This excludes `evergreen` (the single strongest MAT predictor in early runs) because it cannot be determined from a fossil, along with petiole traits and raw tooth counts. Training on traits that can't be measured on fossils would silently undermine the fossil application — the trait set must be identical for training and prediction. If adding predictors, check Dana's list first; this is a hard design requirement, not a preference.

### Tooth trait filling
For confirmed untoothed leaves (cell blank **AND** `margin.score == 1`), tooth trait NAs are biologically real zeros — set them explicitly **before** aggregating to species means in `00_data_cleaning.R`. Without this, tooth traits appear as ~67% NA and get excluded from all models. Tooth-count and tooth-area traits are set to 0; `perim.ratio` is set to 1 (untoothed leaves have a smooth perimeter). Both conditions must be met. Dana requested this for when the model is applied to fossil data.

### Species-level analysis
Modelling is based on species-level trait means (not site means), because phylogeny operates at the species level. Site-level predictions are derived afterward by aggregating species-level predictions within a site.

### Climate targets
MAT and log(MAP) are the primary targets. Other climate variables (coldest month temperature, growing degree days, etc.) covary with MAT/MAP and are not modelled separately.

### Fossil placement
Fossils are commonly known only to family or order level, or may belong to extinct genera. Placement uses a genus → family → order MRCA fallback, with root as the final fallback. Fossil age sets the edge length so the tip sits at the correct time depth. When a fossil predates the crown age of its placement clade, the `PLACEMENT_FALLBACK` flag in `04_fossil_predictions.R` controls behaviour: `"ancestral_branch"` (default, preferred) walks up the tree to the branch alive at the fossil's age and splits it there; `"node"` attaches at the MRCA with a minimal edge.

### Informal fossil taxonomy
Dana's April 2026 fossil dataset marks informal order, family, and genus names with quotation marks. **Never use a quoted rank as phylogenetic evidence in the primary analysis.** `00c_fossil_data_cleaning.R` preserves the reported strings and rank-level informal flags while setting quoted placement ranks to `unknown`. `04_fossil_predictions.R` runs both `formal_only` (primary) and `include_informal` (sensitivity) scenarios and writes placement logs plus a site-level sensitivity table. The legacy output filenames always refer to the formal-only scenario.

### `internal.perimeter.cm` excluded
This raw tooth-linked measurement (~67% NA) is not on Dana's fossil list and remains excluded. The ratio forms (`teeth.perimeter.percm`, `teeth.interior.percm`) are included instead — more appropriate for cross-specimen comparison.

### MAP units
The training data `map` column is in **centimetres**, not millimetres (range ~19–680 cm, mean ~215 cm). All model MAP outputs are therefore in cm, directly comparable to the published DiLP MAP estimates in Peppe et al. (2011). Do not convert.

### Fossil site ages
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

### Species-level vs site-level aggregation — key finding
`04_fossil_predictions.R` compares LM and PIP at both species and site level. Under 10-fold site-grouped CV, **PIP at site level is the best performer for both MAT and MAP** — it takes the top two ranks for each target. PGLS alone (without the PIP covariance correction) is substantially worse — the correction term drives the improvement. PIP acts as a regularization procedure, borrowing signal from phylogenetically close extant relatives. (An earlier claim that "LM site outperforms PIP" came from comparing against published DiLP estimates rather than a CV benchmark; the 10-fold site-grouped CV is the correct evaluation.)

**The recommended variant is target-dependent**: use **impute for MAT** and **complete-case for log(MAP)**. Measured RMSE (full run, 2026-08-06, 1740 species, 93 sites):

| Target | PIP site impute | PIP site cc | best LM site (untoothed-excl, impute) | PGLS alone (best) |
|---|---|---|---|---|
| MAT (°C) | **3.414** | 3.707 | 3.722 | 5.854 |
| log(MAP) | 0.525 | **0.508** | 0.583 | 0.628 |

Note that `specimen_*` and `sp_zero_*` site configs are identical by construction (see issue #10) and report identical RMSE; they are one model, not two.

## Conventions

- Do not add Claude as a co-author in commit messages.
- The `old/` directory contains prior implementations and should not be modified.
- `pw2.a.ratio` is in the fossil trait set but Dana noted "most folks don't measure it" — fossil predictions lacking it rely on imputation.

## R Package Dependencies

`ape`, `phytools`, `caret`, `glmnet`, `ranger`, `pdp`, `gridExtra`, `ggplot2`, `mvtnorm`, `MASS`, `rmarkdown`, `kableExtra`, `dplyr`, `tibble`, `tidyverse` (LMA pipeline), plus the pinned `dilp` fork described above.

`code/Phylogenetically-Informed_Predictions_Source.R` and its self-validation counterpart also carry `library(caper)` / `library(geiger)` calls inherited from the original Freckleton (2015) source; neither package is actually used in that file (no `caper::`, `geiger::`, `comparative.data()`, `fitContinuous()`, or `treedata()` calls), so both `library()` calls are commented out rather than documented as required dependencies.

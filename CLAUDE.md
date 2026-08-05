# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview and Primary Goal

This project predicts paleoclimate (MAT = mean annual temperature, MAP = mean annual precipitation) from fossil leaf morphology using phylogenetically-informed prediction (PIP). **The primary deliverable is fossil climate prediction** — not the extant species model comparison. Every modelling decision should be evaluated against the question: "can this be applied to a fossil leaf?"

The extant species pipeline (scripts 00–03) is a calibration and validation benchmark. Script 04 is the actual scientific product.

## Running the Analysis

Scripts must be run in order. Each sets `setwd()` to a hardcoded path — update line 1 when working on a new machine.

```r
Rscript code/00_data_cleaning.R      # Fill tooth traits, aggregate, build scaffold tree
Rscript code/00c_fossil_data_cleaning.R # Build fossil traits from April 2026 leaf data
Rscript code/01_nophy_regression.R   # Fit non-phylogenetic models (species + site level)
Rscript code/02_phy_regression.R     # Fit PGLS, save pip_components.rds
Rscript code/03_loso_cv.R            # Leave-one-site-out cross-validation
Rscript code/04_fossil_predictions.R # Formal-only primary + informal sensitivity
Rscript code/05_visualizations.R     # Figures from LOSO CV outputs
Rscript code/06_update_dilp_package.R # Regenerate dilp package model objects + constants
```

`code/Phylogenetically-Informed_Predictions_Source.R` is the PGLS library from Freckleton (2015), modified by Gardner et al. (2024).

A parallel `*b*` / LMA pipeline (`00b`–`03b`, `04b_lma_*`, `05b_lma_*`) mirrors these steps to predict leaf mass per area instead of climate; it loads the `dilp` package via `devtools::load_all("~/dilp")` and uses its own calibration inputs.

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
`04_fossil_predictions.R` compares LM and PIP at both species and site level. Under leave-one-site-out CV, **PIP at site level is the best performer for both MAT and MAP** and is the recommended model (impute variant) for fossil climate reconstruction. PGLS alone (without the PIP covariance correction) is substantially worse — the correction term drives the improvement. PIP acts as a regularization procedure, borrowing signal from phylogenetically close extant relatives. (An earlier claim that "LM site outperforms PIP" came from comparing against published DiLP estimates rather than a CV benchmark; the LOSO CV is the correct evaluation.)

## Conventions

- Do not add Claude as a co-author in commit messages.
- The `old/` directory contains prior implementations and should not be modified.
- `pw2.a.ratio` is in the fossil trait set but Dana noted "most folks don't measure it" — fossil predictions lacking it rely on imputation.

## R Package Dependencies

`ape`, `phytools`, `caret`, `glmnet`, `ranger`, `pdp`, `gridExtra`, `ggplot2`, `mvtnorm`, `MASS`, `rmarkdown`, `kableExtra`

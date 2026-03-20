# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is an R-based phylogenetic analysis project studying the relationship between leaf physiognomy (morphological traits) and climate variables (MAT = mean annual temperature, MAP = mean annual precipitation). It compares non-phylogenetic machine learning models against phylogenetically-informed prediction (PIP) models.

## Running the Analysis

Scripts must be run in order. Each script sets `setwd()` to a hardcoded path — update this if running on a different machine.

```r
# Step 1: Clean and prepare species-level data, prune/extend the phylogeny
Rscript code/00_data_cleaning.R

# Step 2: Fit non-phylogenetic models (LM, ElasticNet, RandomForest) via caret
Rscript code/01_nophy_regression.R

# Step 3: Fit phylogenetically-informed PGLS models and run LOOCV
Rscript code/02_phy_regression.R
```

## Pipeline Architecture

**`00_data_cleaning.R`** — Reads the raw Royer leaf morphology/climate CSV and the `best_wcvp.tre_dated` phylogeny. Aggregates specimens to species-level means, then grafts missing species onto the tree at the MRCA of their genus/family/order. Outputs `data/data_species.csv` and `data/tre_pruned.tre`.

**`01_nophy_regression.R`** — Loads cleaned data and fits three model types (LM, ElasticNet via `glmnet`, Random Forest via `ranger`) for MAT and log(MAP) using `caret` with 10-fold CV. Saves fitted models to `models/nophy_models.rds`, exports performance summaries to `tables/`, and PDP plots to `plots/`.

**`02_phy_regression.R`** — Loads `models/nophy_models.rds` and extracts non-zero ElasticNet coefficients to select active traits. Fits a PGLS via `pglmEstLambda()` (ML lambda estimation), then runs a manual leave-one-out cross-validation loop implementing the PIP (Phylogenetically-Informed Predictions) method: for each held-out species, the prediction is the GLS fitted value plus a phylogenetic correction term derived from the VCV matrix partition.

**`code/Phylogenetically-Informed_Predictions_Source.R`** — Core PGLS library sourced by `02_phy_regression.R`. Originally from Freckleton (2015), modified by Gardner et al. (2024). Key functions:
- `pglm()` — fits GLS correcting for phylogeny at a fixed lambda
- `pglmEstLambda()` — ML optimization of Pagel's lambda, then fits `pglm`
- `pglmPredictMissing()` — predicts missing tip values using phylogenetic weights
- `weights.p()` — computes phylogenetic interpolation weights from VCV partitioning

## Key Data Files

| File | Description |
|------|-------------|
| `data/RoyerLeafShapeClimateDataFixedNames_June2012.csv` | Raw leaf morphology + climate data (Royer dataset) |
| `data/best_wcvp.tre_dated` | Dated angiosperm phylogeny (WCVP-matched tip labels: `order_family_genus_species`) |
| `data/data_species.csv` | Species-level averaged traits, output of `00_data_cleaning.R` |
| `data/tre_pruned.tre` | Pruned/extended phylogeny matching species in data |
| `models/nophy_models.rds` | Saved `caret` model objects from `01_nophy_regression.R` |

## R Package Dependencies

`ape`, `caret`, `glmnet`, `ranger`, `pdp`, `gridExtra`, `ggplot2`, `caper`, `geiger`, `mvtnorm`, `MASS`

## Git Commits

Do not add Claude as a co-author in commit messages.

## Key Details

- The phylogeny tip labels use the format `order_family_genus_species`. `00_data_cleaning.R` parses these to build a lookup table for taxonomic grafting.
- Missing species are grafted to the tree by binding at the MRCA of congeners, then confamilials, then conordinals, with edge length set to reach the present.
- The `setwd()` calls use a hardcoded Dropbox path — update line 1 of each script when working on a new machine.
- `data/dat_site.csv` contains site-level averages used for context but not in the main model pipeline.
- The `old/` directory contains prior implementations (Python ML scripts, legacy R PGLS code from collaborators) and should not be modified.

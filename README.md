# phylo_leaf_physiognomy

Phylogenetically-informed analysis of leaf physiognomy as a paleoclimate proxy. This project evaluates whether incorporating evolutionary history improves predictions of mean annual temperature (MAT) and mean annual precipitation (MAP) from leaf morphological traits, benchmarking a phylogenetically-informed prediction (PIP) model against standard machine learning approaches.

## Background

Leaf physiognomy — the relationship between leaf shape/size and climate — is a well-established paleoclimate tool (e.g., DiLP, CLAMP). This project tests whether correcting for shared evolutionary history among species improves prediction accuracy, using the Royer et al. leaf morphology dataset matched to a dated angiosperm phylogeny.

## Pipeline

Run scripts in order. Each script sets `setwd()` to a hardcoded path — update line 1 of each file for your machine.

```r
Rscript code/00_data_cleaning.R      # Aggregate to species, graft onto phylogeny
Rscript code/01_nophy_regression.R   # Fit LM, ElasticNet, Random Forest via caret
Rscript code/02_phy_regression.R     # Fit PGLS + run PIP leave-one-out CV
```

### What each script does

**`00_data_cleaning.R`** — Loads the raw Royer CSV and the WCVP-dated phylogeny (`best_wcvp.tre_dated`). Averages specimens to species level, then grafts species missing from the tree at the MRCA of their genus/family/order. Writes `data/data_species.csv` and `data/tre_pruned.tre`.

**`01_nophy_regression.R`** — Fits three non-phylogenetic models for MAT and log(MAP) using 10-fold CV via `caret`: linear regression (`lm`), elastic net (`glmnet`), and random forest (`ranger`). Saves fitted models to `models/nophy_models.rds`, exports cross-validation metrics to `tables/`, and partial dependence plots to `plots/`.

**`02_phy_regression.R`** — Uses non-zero ElasticNet coefficients from step 2 to select traits, fits a PGLS with ML-estimated Pagel's lambda (`pglmEstLambda`), then runs a manual leave-one-out CV implementing the PIP method: each held-out species is predicted from GLS fitted values plus a phylogenetic correction derived from VCV matrix partitioning. Compares RMSE across all four model types.

## Data

| File | Description |
|------|-------------|
| `data/RoyerLeafShapeClimateDataFixedNames_June2012.csv` | Raw leaf morphology + climate data |
| `data/best_wcvp.tre_dated` | Dated angiosperm phylogeny (tip format: `order_family_genus_species`) |
| `data/smith_and_brown.tre.txt` | Smith & Brown (2018) vascular plant supertree (reference) |

Derived outputs (`data/data_species.csv`, `data/tre_pruned.tre`, `models/`, `plots/`, `tables/`) are excluded from version control — regenerate by running the pipeline.

## Dependencies

```r
install.packages(c("ape", "caret", "glmnet", "ranger", "pdp",
                   "gridExtra", "ggplot2", "caper", "geiger",
                   "mvtnorm", "MASS"))
```

## Key Source File

`code/Phylogenetically-Informed_Predictions_Source.R` — PGLS library originally from Freckleton (2015), modified by Gardner, Baker, Venditti & Organ (2024). Provides `pglm()`, `pglmEstLambda()`, and `pglmPredictMissing()`. Sourced automatically by `02_phy_regression.R`.

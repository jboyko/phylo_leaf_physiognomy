# phylo_leaf_physiognomy

Phylogenetically-informed analysis of leaf physiognomy as a paleoclimate proxy. This project evaluates whether incorporating evolutionary history improves predictions of mean annual temperature (MAT) and mean annual precipitation (MAP) from leaf morphological traits, benchmarking a phylogenetically-informed prediction (PIP) model against standard machine learning approaches.

## Background

Leaf physiognomy — the relationship between leaf shape/size and climate — is a well-established paleoclimate tool (e.g., DiLP, CLAMP). This project tests whether correcting for shared evolutionary history among species improves prediction accuracy, using the Royer et al. leaf morphology dataset matched to a dated angiosperm phylogeny (WCVP).

**The primary goal is fossil prediction.** All modelling decisions — which traits to include, how to handle missing data, how to build the phylogeny — are made with fossil applicability as the priority. The extant species comparison is a benchmark, not the end product.

## Pipeline

Run scripts in order. Each script sets `setwd()` to a hardcoded path — update line 1 of each file for your machine.

```r
Rscript code/00_data_cleaning.R      # Fill tooth traits, aggregate to species, build scaffold tree
Rscript code/01_nophy_regression.R   # Fit LM, ElasticNet, Random Forest via caret
Rscript code/02_phy_regression.R     # Fit PGLS (MAT + MAP), save PIP components
Rscript code/03_comparison.R         # LOOCV, RMSE tables, variable importance, plots
Rscript code/04_fossil_predictions.R # Predict MAT/MAP for fossil specimens via PIP
```

To render the results report:

```r
rmarkdown::render("report.Rmd")
```

### What each script does

**`00_data_cleaning.R`** — Loads the raw Royer CSV and the WCVP-dated phylogeny. Fills tooth trait NAs with 0 (or 1 for `perim.ratio`) for confirmed untoothed leaves (`margin.score == 1`) before aggregating to species-level means and site-level means. Builds a family-level angiosperm backbone (2 crown tips per family across all 515 WCVP families), then grafts training species onto that small scaffold. Outputs `data/data_species.csv`, `data/dat_site.csv`, `data/tre_pruned.tre`, `data/tre_scaffold.tre`, and `data/name_table_full.csv`.

**`01_nophy_regression.R`** — Fits LM, ElasticNet, and Random Forest for MAT and log(MAP) using 10-fold CV via `caret`. **Restricted to the 12 fossil-measurable traits** identified by Dana Royer (pers. comm.) — see trait list below. Pre-imputes once via `bagImpute` before CV for speed. Fits models at both species level and site level. Saves `models/nophy_models.rds` and `models/site_models.rds`.

**`02_phy_regression.R`** — Selects active traits from non-zero ElasticNet coefficients, fits PGLS for MAT and log(MAP) via `pglmEstLambda()`, saves all components for fossil prediction to `models/pip_components.rds`.

**`03_comparison.R`** — All reporting: nophy CV summaries, RF variable importance, PDP plots, vectorised PIP LOOCV (`ŷ₋ᵢ = yᵢ − (Kε)ᵢ/Kᵢᵢ`), RMSE comparison across all models (species-level, site-level, phylogenetic), scatter plots.

**`04_fossil_predictions.R`** — Grafts fossils onto `tre_scaffold.tre`, prunes to training + fossil tips, and applies the PIP formula: `ŷ = Xβ + V[training,fossil] · V[training,training]⁻¹ · ε`. Writes `tables/fossil_predictions.csv`. The `PLACEMENT_FALLBACK` flag at the top of the script controls behaviour when a fossil predates its placement node: `"ancestral_branch"` (default — walks up to split the correct branch at the fossil's age) or `"node"` (attaches at the MRCA with a minimal edge).

**`report.Rmd`** — R Markdown report covering methodology and results. Renders to `report.html`.

## Fossil-Measurable Trait Set

All models are restricted to traits that can be measured from fossil leaf specimens (Dana Royer, pers. comm.). This ensures the extant training model can be directly applied to fossils without imputing unmeasurable characters.

| Trait | Notes |
|-------|-------|
| `ln.leaf.area.mm2` | Log leaf area; Dana's `leaf.area.cm2` (same measurement) |
| `feret.diam.ratio` | Ratio of perpendicular diameters |
| `pw2.a.ratio` | Petiole width² / area; not always measured |
| `margin.score` | 0 = toothed, 1 = untoothed, 0.5 = mixed |
| `perim.ratio` | Set to 1 for untoothed leaves |
| `teeth.perimeter.percm` | Set to 0 for untoothed leaves |
| `teeth.interior.percm` | Set to 0 for untoothed leaves |
| `avt.tooth.area` | Set to 0 for untoothed leaves |
| `tooth.area.blade.area.ratio` | Set to 0 for untoothed leaves |
| `tooth.area.perimeter` | Set to 0 for untoothed leaves |
| `tooth.area.interior` | Set to 0 for untoothed leaves |
| `teeth.blade.area.ratio` | Set to 0 for untoothed leaves |

Traits excluded despite predictive power: `evergreen` (phenological — not determinable from fossils), `major.length.cm`, `petiole.width.cm`, `petiole.area.cm2`, `blade.area.cm2`, `perimeter.cm`, `shape.factor`, `compact`, `perim.area.cm2`, raw tooth counts.

## Data

| File | Description |
|------|-------------|
| `data/RoyerLeafShapeClimateDataFixedNames_June2012.csv` | Raw leaf morphology + climate data |
| `data/best_wcvp.tre_dated` | Dated angiosperm phylogeny (tip format: `order_family_genus_species`) |
| `data/fossil_traits.csv` | User-provided fossil specimen traits for prediction |

Derived outputs (`data/data_species.csv`, `data/tre_pruned.tre`, `data/tre_scaffold.tre`, `data/name_table_full.csv`, `models/`, `plots/`, `tables/`) are excluded from version control — regenerate by running the pipeline.

## Dependencies

```r
install.packages(c("ape", "phytools", "caret", "glmnet", "ranger",
                   "pdp", "gridExtra", "ggplot2", "mvtnorm", "MASS",
                   "rmarkdown", "kableExtra"))
```

## Key Source File

`code/Phylogenetically-Informed_Predictions_Source.R` — PGLS library originally from Freckleton (2015), modified by Gardner, Baker, Venditti & Organ (2024). Provides `pglm()`, `pglmEstLambda()`, and `pglmPredictMissing()`. Sourced automatically by `02_phy_regression.R`.

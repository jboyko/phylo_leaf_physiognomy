# phylo_leaf_physiognomy

Phylogenetically-informed analysis of leaf physiognomy as a paleoclimate proxy. This project evaluates whether incorporating evolutionary history improves predictions of mean annual temperature (MAT) and mean annual precipitation (MAP) from leaf morphological traits, benchmarking a phylogenetically-informed prediction (PIP) model against non-phylogenetic linear regression (LM), then applies the fitted models to fossil leaf floras.

## Background

Leaf physiognomy — the relationship between leaf shape/size and climate — is a well-established paleoclimate tool (e.g., DiLP, CLAMP). This project tests whether correcting for shared evolutionary history among species improves prediction accuracy, using the Royer et al. leaf morphology dataset matched to a dated angiosperm phylogeny (WCVP).

**The primary goal is fossil prediction.** All modelling decisions — which traits to include, how to handle missing data, how to build the phylogeny — are made with fossil applicability as the priority. The extant species comparison is a benchmark, not the end product.

## Pipeline

Run scripts in order from the repository root. Each script sources `code/setup.R` as its first line, which locates the repository root automatically and creates the `models/`, `tables/`, `plots/` output directories — no per-machine path edits needed.

```r
Rscript code/00_data_cleaning.R      # Fill tooth traits, aggregate to species, build scaffold tree
Rscript code/00c_fossil_data_cleaning.R # Build fossil traits from the April 2026 leaf-level data
Rscript code/01_nophy_regression.R   # Fit LM via caret (species + site level)
Rscript code/02_phy_regression.R     # Fit PGLS (MAT + MAP), save PIP components
Rscript code/03_loso_cv.R            # 10-fold site-grouped CV across all 12 model configs, RMSE tables, model summaries
Rscript code/04_fossil_predictions.R # Predict MAT/MAP under both taxonomy scenarios
Rscript code/05_visualizations.R     # Figures from the 10-fold site-grouped CV outputs
Rscript code/06_update_dilp_package.R # Maintenance utility: regenerate dilp package model objects + constants
Rscript code/07_type2_degradation.R  # Sensitivity: adjustment degradation vs. placement distance
Rscript code/08_type1_degradation.R  # Sensitivity: error from coarser (genus/family/order) placement
Rscript code/09_phylo_adj_field.R    # Map the phylogenetic adjustment field over tre_pruned
```

A parallel `*b*` / LMA pipeline (`00b`–`03b`, `04b_lma_*`, `05b_lma_*`) mirrors
`00`–`05` to predict leaf mass per area instead of climate. It loads the `dilp`
package via `devtools::load_all("~/dilp")` and uses its own calibration inputs
(`data/extra_calibration_data_for_LMA.csv`, `data/lma_species.csv`).

### What each script does

**`00_data_cleaning.R`** — Loads the raw Royer CSV and the WCVP-dated phylogeny. Fills tooth trait NAs with 0 (or 1 for `perim.ratio`) for confirmed untoothed leaves (`margin.score == 1`) before aggregating. Outputs three site-level datasets for model comparison: `dat_site.csv` and `dat_site_sp_zero.csv` (both morphotypes → site, zero-filled, identical in the current implementation), and `dat_site_untoothed_excl.csv` (morphotypes → site, untoothed excluded from tooth-trait averages). Builds a family-level angiosperm backbone (2 crown tips per family across all 515 WCVP families), then grafts training species onto that small scaffold. A `BUILD_PHYLOGENY` flag (default `TRUE`) skips the slow tree-building sections when set to `FALSE` and only data changes. Outputs `data/data_species.csv`, the three site CSVs, `data/tre_pruned.tre`, `data/tre_scaffold.tre`, and `data/name_table_full.csv`.

**`00c_fossil_data_cleaning.R`** — Converts Dana Royer's April 2026 leaf-level fossil dataset into one row per species × site using `dilp`. It restores the Palacio de los Loros PL1/PL2 assignments, applies the documented site ages, preserves the reported taxonomy and quote flags, and writes formal-only placement columns to `data/fossil_traits.csv`. Quoted genus/family/order names are informal: the primary analysis censors the quoted rank, while the sensitivity analysis can use its unquoted value provisionally.

**`01_nophy_regression.R`** — Fits **LM only** (via `caret`) at **species level** for MAT and log(MAP); earlier ElasticNet and Random Forest comparisons are not part of the current analysis. **Restricted to the 12 fossil-measurable traits** identified by Dana Royer (pers. comm.) — see trait list below. Also fits LM at site level across six combinations (3 datasets × bagImpute / complete-case). All site configs stored under `site_models$configs`; backward-compatible top-level keys preserved. Saves `models/nophy_models.rds` and `models/site_models.rds`.

**`02_phy_regression.R`** — Uses the same 12 predictor set as the LM (loaded from `nophy_models.rds`) — no ENet-based trait selection, ensuring fair comparison. Fits PGLS for MAT and log(MAP) via `pglmEstLambda()` in two variants: `impute` (bagImpute) and `cc` (complete-case, phylogeny pruned to complete species). All variants stored under `pip_components$configs`; backward-compatible top-level keys point to the `impute` variant. Saves `models/pip_components.rds`.

Species-level training rows are currently operational labels: named species are averaged across sites, while records lacking a species epithet are pooled under genus-level keys. Those genus-only pools can combine several sites or morphotypes and await taxonomic review. They are retained for the current analysis; changing this convention requires rebuilding the training data and rerunning full CV.

**`04_fossil_predictions.R`** — Grafts each fossil occurrence onto `tre_scaffold.tre` at its reported site age and runs four prediction approaches under two taxonomy scenarios. The canonical `tables/fossil_site_comparison.csv` and `tables/fossil_predictions.csv` use formal taxonomy only: quoted genus, family, and order ranks are excluded from placement evidence. Scenario-specific predictions, placement logs, and `tables/fossil_taxonomy_sensitivity.csv` compare that primary interpretation with provisional inclusion of quoted taxonomy:

| Method | Training data | Fossil input | Aggregation |
| ------ | ------------ | ------------ | ----------- |
| LM sp | extant species means | fossil species grand means | average per site |
| LM site | extant site means | fossil site-mean traits | direct (one per site) |
| PIP sp | extant species (PGLS) | fossil species grand means | average per site |
| PIP+site | extant species (PGLS) | site-specific species means | average per site |

PIP sp uses fossil species grand means; PIP+site uses site-specific species means and occurrence-specific placements. The `PLACEMENT_FALLBACK` flag controls behaviour when a fossil predates its placement node: `"ancestral_branch"` (default) or `"node"`.

**`03_loso_cv.R`** — 10-fold cross-validation, grouped by site, across all 12 model configurations. Sites are ranked by site MAT and assigned round-robin to folds (roughly 9 sites held out per fold); each fold retrains all models from scratch on the remaining sites and predicts the held-out sites. For species-level MAP models, log predictions are averaged within site for evaluation; exponentiating that mean gives the geometric site MAP estimate, matching the fossil MAP estimator. Whole sites are withheld, though their extant species can occur at training sites; the benchmark therefore tests site transfer with known extant relationships, not fossil-placement accuracy. It is **not** leave-one-site-out (that would be 93 folds, one per site) — see the naming note below. Saves per-fold model objects to `models/loso_cv_fold_XX.rds`. Outputs: `tables/loso_cv_rmse.csv`, `tables/loso_cv_site_predictions.csv`, `tables/loso_cv_model_coefs.csv` (per-fold coefficient estimates), `tables/loso_cv_model_fit.csv` (per-fold R², residual SE, lambda).

**`05_visualizations.R`** — Builds the figures and diagnostic tables from the `03_loso_cv.R` outputs (RMSE comparison, observed-vs-predicted, coefficient stability, Pagel's lambda across folds).

**`06_update_dilp_package.R`** — Maintenance utility (not part of the main pipeline). Regenerates the pre-baked model objects and RMSE constants shipped inside the `dilp` package from the current CV results. Needs a local checkout of the `dilp` source.

**`07_type2_degradation.R`** and **`08_type1_degradation.R`** — Sensitivity analyses on the PIP adjustment: how much it degrades with placement distance from the true position (07), and how much error is introduced by only knowing a fossil to genus/family/order rather than species (08).

**`09_phylo_adj_field.R`** — Visualizes the PIP phylogenetic-adjustment field over every tip and internal node of `tre_pruned`, rendered as fan phylogenies with order-level bands.

## A note on cross-validation naming

`03_loso_cv.R` and its outputs are named `loso` for continuity, but the procedure
is **10-fold cross-validation with sites as the grouping unit**, not
leave-one-site-out. The 93 calibration sites are ranked by site MAT and assigned round-robin to
10 folds, so roughly 9 sites are held out per fold. Whole sites are always held
out together and every model is refitted from scratch on the remaining sites, so
no specimen from a held-out site informs its own prediction. Figure titles
produced by `05_visualizations.R` have been updated to say "10-fold site-grouped
CV"; the `loso_cv_*` file, column, and object names remain unchanged for
continuity.

## 10-Fold Site-Grouped CV Model Configurations

`03_loso_cv.R` evaluates 6 model types, each in two NA-handling variants (bagImpute and complete-case), for a total of 12 configurations per climate target. Because the held-out unit is always a site, all predictions are site-level values. Column names use `sp_site` to indicate a species-level model whose predictions are averaged to site, and `site` to indicate a model trained directly on site-level data.

**Species-level LM (`lm_sp_site`)**
Ordinary least squares trained on species grand means (averaged across all training sites). Each held-out species gets a prediction from Xβ; those are aggregated to a site estimate. No use of phylogenetic information. Averaging across training sites loses within-species climatic and trait variation, whereas held-out and fossil prediction rows preserve within-site species means; this is a limitation of this calibration implementation.

**Site-level LM — `specimen` (`lm_site_specimen`)**
Ordinary least squares trained on site means computed by averaging `dilp` morphotype means within sites. In the current data this is identical by construction to `sp_zero`; it does not weight raw specimens directly.

**Site-level LM — `sp_zero` (`lm_site_sp_zero`)**
Ordinary least squares trained on the current `dilp` morphotype means within site, with zero-filling for confirmed untoothed leaves. It is identical by construction to the `specimen` configuration in the present data.

**Site-level LM — `untoothed_excl` (`lm_site_untoothed_excl`)**
Same two-step aggregation as `sp_zero`, but no zero-filling — untoothed species contribute `NA` to tooth traits, which are dropped via `na.rm = TRUE`. This local aggregation configuration is distinct from fixed published DiLP coefficients.

**PGLS (`pgls_sp_site`)**
Generalised least squares on species means with the phylogenetic VCV as the covariance structure. Beta is estimated as (X'V⁻¹X)⁻¹X'V⁻¹y, accounting for phylogenetic non-independence during fitting. Prediction for held-out species is Xβ only — no phylogenetic correction is applied at prediction time. Species predictions are averaged to site.

**PIP (`pip_sp_site`)**
Shares the same beta estimation as PGLS. Adds a phylogenetic covariance correction at prediction time: ŷ = Xβ + C·K·ε, where C is the cross-covariance between held-out and training species, K = V⁻¹, and ε are training residuals. Borrows signal from phylogenetically close training species. Species predictions are averaged to site.

**Interpreting model comparisons:** PGLS and PIP share fitted coefficients, so their difference isolates the prediction-time covariance correction. LM species versus PIP also changes coefficient estimation, while site-level LM versus PIP changes fitting and aggregation level. The requested LM site sp+zero comparison is a practical baseline, not a pure test of phylogeny alone.

| Axis | Options |
| --- | --- |
| Training level | species grand mean (`lm_sp_site`, `pgls_sp_site`, `pip_sp_site`), site mean (`lm_site_*`) |
| Site aggregation method | specimen direct, sp_zero, peppe (site LM only) |
| Phylogenetic adjustment | none (LM), fit-only (PGLS), fit + predict (PIP) |
| NA handling | bagImpute, complete-case |

## Fossil-Measurable Trait Set

All models are restricted to traits that can be measured from fossil leaf specimens (Dana Royer, pers. comm.). This ensures the extant training model can be directly applied to fossils without imputing unmeasurable characters.

| Trait | Notes |
| ----- | ----- |
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
| ---- | ----------- |
| `data/Peppe_2011_calibration_data_leaf_level_clean.csv` | Raw leaf morphology + climate data (read by `00_data_cleaning.R`) |
| `data/best_wcvp.tre_dated` | Dated angiosperm phylogeny (tip format: `order_family_genus_species`) |
| `data/Peppe_2011_fossil_data_April_2026_leaf_level_clean.csv` | Dana Royer's updated leaf-level fossil measurements and taxonomy |
| `data/fossil_traits.csv` | Tracked, reproducible fossil analysis input generated by `00c_fossil_data_cleaning.R` — one row per species × site, with formal-only placement taxonomy, reported taxonomy, informal-rank flags, site age, and the 12 trait columns. MAP units throughout are **cm** (same as the training data). |
| `data/extra_calibration_data_for_LMA.csv`, `data/lma_species.csv` | Calibration inputs specific to the `*b*` / LMA pipeline |
| `data/RoyerLeafShapeClimateDataFixedNames_June2012.csv` | Original Royer leaf-shape + climate dataset. **Superseded** by `Peppe_2011_calibration_data_leaf_level_clean.csv` and read by no script; retained as the provenance record for the calibration set (it carries the pre-filtering columns, e.g. `evergreen`, petiole traits, `internal.perimeter.cm`). |

Other derived outputs (`data/data_species.csv`, `data/tre_pruned.tre`, `data/tre_scaffold.tre`, `data/name_table_full.csv`, `models/`, `plots/`, `tables/`) are excluded from version control — regenerate by running the pipeline.

## Dependencies

```r
install.packages(c("ape", "phytools", "caret", "glmnet", "ranger",
                   "pdp", "gridExtra", "ggplot2", "mvtnorm", "MASS",
                   "rmarkdown", "kableExtra", "dplyr", "tidyr", "tibble"))

# Pinned dilp fork — the CRAN release carries different regression
# constants and will not reproduce these results:
remotes::install_github("jboyko/dilp", ref = "b29be909355f7edb7809bf718f189b7a53b789d2")
```

`caret`'s bagged-tree imputation additionally needs the `ipred` and `e1071` packages installed (not loaded directly).

`code/Phylogenetically-Informed_Predictions_Source.R` carries `library(caper)` / `library(geiger)` calls inherited from the original Freckleton (2015) source, but neither package is used in that file — both calls are commented out, so `caper` and `geiger` are not required.

## Key Source File

`code/Phylogenetically-Informed_Predictions_Source.R` — PGLS library originally from Freckleton (2015), modified by Gardner, Baker, Venditti & Organ (2024). Provides `pglm()`, `pglmEstLambda()`, and `pglmPredictMissing()`. Sourced automatically by `02_phy_regression.R`.


## September 2026 climate update

Primary fossil predictions now use fixed extant scaffold anchors and retain all 361 occurrences. Tests cover insertion-order invariance, site ages, and numerical agreement with the updated local dilp source. Run `Rscript tests/test_fossil_placement.R`, `Rscript tests/test_site_prediction.R`, and `DILP_SOURCE=/path/to/updated/dilp Rscript tests/test_dilp_parity.R` after generating models and fossil outputs.

The pinned April dilp commit remains the raw `dilp()` preprocessing dependency. It does **not** contain the revised `dilp_pgls()` implementation. The September package changes are local and must be published/pinned separately before distributing a package-based reproduction claim. `code/06_update_dilp_package.R` stages fitted trait imputers, coefficients, covariance weights, lambda, and current validation error scales into `dilp_update/`.

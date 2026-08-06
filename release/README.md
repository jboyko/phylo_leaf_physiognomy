# Phylogenetically-informed prediction of paleoclimate from leaf physiognomy

Reproduction package for the analysis. It predicts mean annual temperature (MAT)
and mean annual precipitation (MAP) from leaf morphology using phylogenetically
informed prediction (PIP), and applies the fitted models to the Peppe et al.
(2011) fossil floras. A parallel pipeline predicts leaf mass per area (LMA)
instead of climate.

Every model is restricted to the twelve traits that can be measured on a fossil
leaf, so the trait set used for training is identical to the one used for
prediction.

## Quick start

```bash
cd release
Rscript setup.R      # install dependencies (once)
Rscript run_all.R    # run the whole pipeline
```

Scripts are also runnable individually and in any order that respects the
dependency chain below. All of them must be run **from the release root**; they
use relative paths and will stop with an error otherwise.

Expect several hours end to end on a laptop. Measured on an Apple M-series
machine: `00_data_cleaning.R` ~4 min (grafting species onto the angiosperm
backbone), `01_nophy_regression.R` ~3 min, `02_phy_regression.R` ~5 min,
`04_fossil_predictions.R` ~1 min. `03_loso_cv.R` dominates the total, since it
refits every model from scratch in each of ten folds.

`models/` grows to a few hundred MB — the fitted `caret` objects and the per-fold
CV objects are large. Allow ~2 GB of free disk space.

## Dependencies

`setup.R` installs everything. The pipeline needs R ≥ 4.1 (it uses the native
`|>` pipe and `\(x)` lambda syntax) plus:

| Source | Packages |
| --- | --- |
| CRAN | `ape`, `phytools`, `caper`, `geiger`, `caret` (with `ipred`, `e1071` for `bagImpute`), `dplyr`, `tidyr`, `tibble`, `ggplot2`, `gridExtra`, `mvtnorm`, `MASS` |
| GitHub | `jboyko/dilp`, pinned to commit `b29be90` |

`dilp` is pinned deliberately. The CRAN release of the package carries different
regression constants and will not reproduce these results. `setup.R` installs the
pinned commit; `code/_setup.R` refuses to run if it is missing.

## Pipeline

`code/_setup.R` is sourced by every script. It verifies the working directory,
creates `models/`, `tables/` and `plots/`, and sets the RNG seed to 42 so the
bagged-tree imputation and the CV fold assignment are reproducible.

### Climate pipeline

| Script | Reads | Writes |
| --- | --- | --- |
| `00_data_cleaning.R` | `Peppe_2011_calibration_data_leaf_level_clean.csv`, `best_wcvp.tre_dated` | `data/data_species.csv`, `data/dat_site*.csv`, `data/tre_scaffold.tre`, `data/tre_pruned.tre`, `data/name_table_full.csv` |
| `00c_fossil_data_cleaning.R` | `Peppe_2011_fossil_data_April_2026_leaf_level_clean.csv` | `data/fossil_traits.csv` |
| `01_nophy_regression.R` | `data/data_species.csv`, `data/dat_site*.csv` | `models/nophy_models.rds`, `models/site_models.rds` |
| `02_phy_regression.R` | `data/tre_pruned.tre`, `models/nophy_models.rds` | `models/pip_components.rds` |
| `03_loso_cv.R` | raw calibration data, `data/tre_pruned.tre` | `tables/loso_cv_*.csv`, `models/loso_cv_fold_XX.rds` |
| `04_fossil_predictions.R` | `models/*.rds`, `data/tre_scaffold.tre`, `data/fossil_traits.csv` | `tables/fossil_*.csv` |
| `05_visualizations.R` | `tables/loso_cv_*.csv` | `plots/fig1–4`, `tables/prediction_diagnostics.csv` |

**`00_data_cleaning.R`** processes specimen data through `dilp()`, fills tooth
traits for confirmed untoothed leaves (zero, or 1 for `perim.ratio`) before
aggregation, and aggregates to species means and to three site-level variants.
It then builds a family-level angiosperm backbone from the dated WCVP tree — two
crown tips per family — and grafts the training species onto it. `BUILD_PHYLOGENY
<- FALSE` skips the tree-building sections when only the trait data has changed.

**`00c_fossil_data_cleaning.R`** builds one row per fossil species × site,
restores the Palacio de los Loros PL1/PL2 split, applies the site ages in
`code/fossil_taxonomy.R`, and records which taxonomic ranks were reported
informally (in quotation marks).

**`01_nophy_regression.R`** fits ordinary least squares at species level and at
site level, in bagged-imputation and complete-case variants. Only LM is fitted —
earlier elastic-net and random-forest comparisons are not part of this analysis.

**`02_phy_regression.R`** fits PGLS with `pglmEstLambda()` on the same predictor
set, and stores the components PIP needs at prediction time: β, λ, the
λ-transformed VCV, and the training residuals.

**`03_loso_cv.R`** is the predictive benchmark. See the naming note below.

**`04_fossil_predictions.R`** grafts each fossil species onto the scaffold tree
at its genus, family or order MRCA, using the fossil's age to set the tip depth,
and predicts under four methods:

| Method | Training data | Fossil input |
| --- | --- | --- |
| LM sp | extant species means | fossil species grand means, averaged to site |
| LM site | extant site means | fossil site-mean traits, one prediction per site |
| PIP sp | extant species (PGLS) | fossil species grand means, averaged to site |
| PIP+site | extant species (PGLS) | site-specific species means, averaged to site |

It runs two taxonomy scenarios. `formal_only` is primary and censors any rank
reported in quotation marks; `include_informal` is a sensitivity analysis that
takes the quoted name at face value. The unsuffixed output files
(`tables/fossil_predictions.csv`, `tables/fossil_site_comparison.csv`) always
hold the formal-only results.

`PLACEMENT_FALLBACK` controls what happens when a fossil predates the crown age
of its placement clade. The default `"ancestral_branch"` walks up the tree to the
branch alive at the fossil's age and splits it there; `"node"` attaches at the
MRCA with a minimal edge instead.

### Degradation and adjustment-field analyses

| Script | Question |
| --- | --- |
| `07_type2_degradation.R` | How much does the PIP adjustment degrade as the placement moves further from the true position? Analytic leave-one-out over all species pairs. |
| `08_type1_degradation.R` | How much error comes from knowing a fossil only to genus, family or order rather than to species? |
| `09_phylo_adj_field.R` | Maps the phylogenetic adjustment a hypothetical fossil would receive at every node of the training tree. |

### LMA pipeline

`00b`–`05b` mirror the climate pipeline but predict `log10(LMA)` from petiole
width² / area and the same fossil-measurable trait set. They add
`extra_calibration_data_for_LMA.csv` to the calibration data and build their own
pruned tree. `00b` reuses `data/tre_scaffold.tre`, so `00_data_cleaning.R` must
be run first.

### Not part of the analysis

`06_update_dilp_package.R` regenerates the pre-baked model objects and constants
that ship inside the `dilp` package. It is a maintenance utility, needs a local
checkout of the `dilp` source, and is excluded from `run_all.R`. Point it at a
checkout with `DILP_SYSDATA=/path/to/dilp/R/sysdata.rda`.

## A note on cross-validation naming

`03_loso_cv.R` and its outputs are named `loso` for continuity, but the procedure
is **10-fold cross-validation with sites as the grouping unit**, not
leave-one-site-out. The 93 calibration sites are ranked by site MAT and dealt into
10 folds, so roughly 9 sites are held out per fold. Whole sites are always held
out together and every model is refitted from scratch on the remaining sites, so
no specimen from a held-out site informs its own prediction. Figure titles
produced by `05_visualizations.R` carry the same `LOSO` label.

## Model configurations

`03_loso_cv.R` evaluates six model types, each in a bagged-imputation and a
complete-case variant. Because the held-out unit is a site, all predictions are
site-level values.

| Key | Description |
| --- | --- |
| `lm_sp_site` | OLS on species grand means; species predictions averaged to site |
| `lm_site_specimen` | OLS on site means from direct specimen averaging |
| `lm_site_sp_zero` | OLS on site means via specimens → species-within-site → site |
| `lm_site_peppe` | As `sp_zero` but without tooth zero-filling, matching Peppe et al. (2011) |
| `pgls_sp_site` | GLS on species means with the phylogenetic VCV; prediction is Xβ only |
| `pip_sp_site` | Same fit as PGLS, plus the prediction-time correction ŷ = Xβ + C·K·ε |

Column names follow
`{method}_{train_level}_{pred_level}_{config}_{na_handling}_{target}`, e.g.
`pip_sp_site_impute_mat`.

PGLS and PIP share identical coefficients; they differ only in whether the
covariance correction is applied at prediction time. Any performance difference
between them is therefore attributable to that correction alone.

## Fossil-measurable trait set

| Trait | Note |
| --- | --- |
| `ln.leaf.area.mm2` | Log leaf area |
| `feret.diam.ratio` | Ratio of perpendicular diameters |
| `pw2.a.ratio` | Petiole width² / area; not measured by all workers, so often imputed for fossils |
| `margin.score` | 0 = toothed, 1 = untoothed, 0.5 = mixed |
| `perim.ratio` | Set to 1 for untoothed leaves |
| `teeth.perimeter.percm` | Set to 0 for untoothed leaves |
| `teeth.interior.percm` | Set to 0 for untoothed leaves |
| `avt.tooth.area` | Set to 0 for untoothed leaves |
| `tooth.area.blade.area.ratio` | Set to 0 for untoothed leaves |
| `tooth.area.perimeter` | Set to 0 for untoothed leaves |
| `tooth.area.interior` | Set to 0 for untoothed leaves |
| `teeth.blade.area.ratio` | Set to 0 for untoothed leaves |

Traits excluded because they cannot be scored on a fossil, despite predictive
power on extant material: `evergreen` (the strongest single MAT predictor),
petiole dimensions, and raw tooth counts. A trait absent from a fossil for a
different reason — because the worker did not measure it — is imputed with the
bagged-tree model fitted on the extant training data.

For confirmed untoothed leaves (`margin.score == 1`) the missing tooth traits are
biologically real zeros rather than missing data, and are set explicitly before
species aggregation. Without this they appear as roughly 67% NA and are dropped
from every model.

## Units

The `map` column is in **centimetres**, not millimetres (range ~19–680 cm). All
model MAP output is therefore in cm, directly comparable to the published DiLP
estimates in Peppe et al. (2011).

## Fossil site ages

Midpoints of the published ranges, used to set tip depth during phylogenetic
placement. Defined in `code/fossil_taxonomy.R`.

| Site | Age (Ma) |
| --- | --- |
| Fox Hills | 66.5 |
| Williston Basin I | 64.75 |
| Williston Basin II | 63.5 |
| Williston Basin III | 59.75 |
| Palacio de los Loros PL1 and PL2 | 61.7 |
| Cerrejón | 58.0 |
| Hubble Bubble | 55.8 |
| Laguna del Hunco | 51.9 |
| Republic | 49.4 |
| Bonanza | 47.3 |

## Data

| File | Description |
| --- | --- |
| `data/Peppe_2011_calibration_data_leaf_level_clean.csv` | Extant leaf-level morphology and climate, 93 sites |
| `data/Peppe_2011_fossil_data_April_2026_leaf_level_clean.csv` | Fossil leaf-level measurements and taxonomy, 11 sites |
| `data/extra_calibration_data_for_LMA.csv` | Additional petiole-width records used only by the LMA pipeline |
| `data/best_wcvp.tre_dated` | Dated angiosperm phylogeny; tip labels are `order_family_genus_species` |

Everything else under `data/` is derived and regenerated by `00_data_cleaning.R`,
`00b_data_cleaning.R` and `00c_fossil_data_cleaning.R`.

## Repository layout

```
release/
  setup.R          install pinned dependencies
  run_all.R        run the pipeline end to end
  code/            analysis scripts and shared helpers
  data/            the four raw input files
  tests/           assertions on the fossil taxonomy rules
  models/          fitted model objects (created on run)
  tables/          CSV results (created on run)
  plots/           figures (created on run)
```

Run the tests with `Rscript tests/test_fossil_taxonomy.R` from the release root.
They check the informal-taxonomy rules, the site ages and the shape of
`data/fossil_traits.csv`, so `00c_fossil_data_cleaning.R` must have run first.

## Third-party code

`code/Phylogenetically-Informed_Predictions_Source.R` is the PGLS library
originally written by R. P. Freckleton (2015), used in Franks et al. (2012), and
modified and annotated by Gardner, Baker, Venditti and Organ (2024). It provides
`pglm()`, `pglmEstLambda()` and `pglmPredictMissing()`. It is redistributed here
unmodified.

## References

Peppe, D. J. et al. (2011) Sensitivity of leaf size and shape to climate:
global patterns and paleoclimatic applications. *New Phytologist* 190, 724–739.

Lowe, A. J. et al. (2024) Digital leaf physiognomy. Methods implemented in the
`dilp` R package.

Gardner, E. E., Baker, J., Venditti, C. and Organ, C. (2024)
Phylogenetically-informed predictions.

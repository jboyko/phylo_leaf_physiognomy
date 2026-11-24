# PART 4: Empirical Application

## Methods

### Fossil leaf data

Leaf morphological measurements were taken from 10 late Cretaceous through early Eocene fossil sites from Peppe et al. (2011), spanning ages from 67 Ma (Fox Hills) to 47.3 Ma (Bonanza). Palacio de los Loros comprises the PL1 and PL2 collections, combined before species means are calculated. The dataset contains 360 species-by-site rows, with one row per fossil species at each site at which it occurs. Each row records the same 12 fossil-measurable predictors used to train the climate models described in Part 1. Each fossil species entry also carries four taxonomic fields (species, genus, family, order) and the site age in millions of years (age_ma). These fields are required for phylogenetic placement and are the same columns expected by the `dilp_pgls()` function in the dilp R package. A subset of fossil species also carry measurements of PW²/A (petiole width squared divided by blade area), the single predictor used in the LMA model; species without this measurement were excluded from LMA prediction.

### Phylogenetic placement

Each fossil occurrence was grafted onto the angiosperm scaffold tree (described in Part 1) using a genus-then-family-then-order fallback hierarchy and that occurrence's site age. The primary analysis skips any rank marked with quotation marks or a trailing question mark. The sensitivity analysis tries those ranks provisionally after removing the markers. If two or more extant scaffold tips match the highest available formal rank, their MRCA is the placement target; a single matching tip defines the target branch; otherwise the next formal rank is used, with the root as the final fallback. The edge length is chosen so that the fossil tip falls at its correct time depth. When a fossil predates the crown age of its placement clade, the code walks up the tree to find a branch alive at the fossil's age and splits that branch. This ancestral-branch approach ensures that no fossil is placed at a node younger than itself.

### PIP prediction

After placement, the phylogenetic variance-covariance matrix was computed for the combined set of extant training species and fossil occurrences using covariance blocks from the grafted tree with the climate scaffold root preserved. The cross-covariance block $V_{cross}$ between the $n$ training species and the $m$ fossil occurrences was extracted and scaled by Pagel's $\lambda$ estimated during PGLS fitting. The phylogenetic adjustment $V_{cross}^T V_{inv} e$ was computed as described in Part 1, with $V_{inv}$ and $e$ taken from the fitted PGLS components stored during model training. The trait-based component $X\beta$ was computed using site-specific trait means for each fossil occurrence. Each occurrence is also placed at its own age, so both its cross-covariance and phylogenetic adjustment can differ when the same label occurs at more than one site or age. Missing trait values were filled by bagged-tree imputation before constructing the design matrix. Predictions on the log(MAP) scale were averaged within each site and then exponentiated, giving geometric site means. Site-level MAT estimates remain arithmetic means of the occurrence-level predictions.

LMA predictions followed the same PIP framework using a separate model fitted to the extant LMA calibration dataset. The single predictor is log₁₀(PW²/A) (Dana Royer pers. comm.). For each fossil species with a valid PW²/A measurement, the design matrix $X$ contains an intercept column and a column of log₁₀(PW²/A) values aggregated to within-site species means. No imputation was applied; species without PW²/A measurements were excluded. Each LMA occurrence is placed at its own site age, with $V_{cross}$, $V_{inv}$, and $e$ drawn from the LMA model’s fitted components. Predictions on the log₁₀(LMA) scale were back-transformed as $10^{\hat{y}}$ after averaging occurrence predictions within each site, yielding geometric means in g m⁻².

The climate predictions can be reproduced using the pinned `dilp_pgls()` implementation. The function accepts a specimen-level data frame with the standard DiLP trait columns plus the species, genus, family, order, and age_ma fields, and returns site-level and species-level predictions together with a placement log recording how each fossil species was grafted onto the scaffold phylogeny.

## Results

All 360 fossil occurrences were retained in the formal-only climate analysis. PIP estimates range from 14.1 to 21.6 °C for MAT and 146 to 216 cm for MAP. The practical comparison uses an imputed twelve-trait site regression.

**Table 6.** Site-level climate estimates for ten fossil floras. Brackets give nominal 95% PIP prediction intervals with a shared site discrepancy term. MAP is back-transformed from the site mean and interval on the log scale. The LM predicts directly at site level. n counts fossil occurrences.

| Site | Age (Ma) | MAT PIP (°C) [interval] | MAT LM (°C) | MAP PIP (cm) [interval] | MAP LM (cm) | n |
| --- | --- | --- | --- | --- | --- | --- |
| Fox Hills | 67 | 16.5 [7.8, 25.1] | 18.1 | 154 [50, 472] | 83 | 24 |
| Williston Basin I | 65.07 | 17.3 [7.1, 27.5] | 14.5 | 169 [51, 561] | 137 | 20 |
| Palacio de los Loros | 64.08 | 17.2 [9.4, 25.0] | 13.5 | 166 [57, 489] | 108 | 29 |
| Williston Basin II | 63.8 | 16.5 [7.4, 25.6] | 13.7 | 162 [52, 507] | 102 | 23 |
| Williston Basin III | 60.4 | 16.5 [8.2, 24.8] | 16.0 | 161 [53, 484] | 107 | 18 |
| Cerrejon | 59 | 21.6 [12.6, 30.7] | 25.9 | 216 [70, 673] | 283 | 45 |
| Hubble Bubble | 55.8 | 19.2 [10.9, 27.4] | 17.8 | 166 [55, 500] | 96 | 16 |
| Laguna del Hunco | 52 | 16.7 [6.3, 27.1] | 11.0 | 168 [50, 558] | 139 | 119 |
| Republic | 51.18 | 14.1 [5.9, 22.3] | 8.6 | 146 [49, 437] | 67 | 41 |
| Bonanza | 47.3 | 17.1 [8.6, 25.6] | 10.3 | 149 [49, 451] | 116 | 25 |

The site-discrepancy intervals covered 96.7% of extant MAT and 93.5% of extant log(MAP) observations in a fold-exclusion check with overlapping training sets. Fossil intervals condition on the fitted model, selected placements, ages, and imputed traits. The primary placements were 195 at root, 105 at family, 37 at genus, and 23 at order. Provisionally including uncertain ranks changed rounded site PIP estimates by up to 0.6 °C and 4 cm.

PIP geometric site means range from 73.1 to 113.6 g m⁻². The analysis includes 192 species occurrences across 9 sites with usable PW²/A. Palacio de los Loros contributes 23 occurrences (Table 7).

**Table 7.** PIP and LM geometric mean LMA in g m⁻² at nine fossil sites. n counts species occurrences with usable PW²/A.

| Site | Age (Ma) | LMA PIP (g m⁻²) | LMA LM (g m⁻²) | n |
| --- | --- | --- | --- | --- |
| Fox Hills | 67 | 87.7 | 85.1 | 16 |
| Williston Basin I | 65.07 | 84.8 | 81.7 | 11 |
| Palacio de los Loros | 64.08 | 83.0 | 80.9 | 23 |
| Williston Basin II | 63.8 | 77.9 | 77.2 | 13 |
| Williston Basin III | 60.4 | 73.1 | 72.2 | 13 |
| Hubble Bubble | 55.8 | 90.9 | 96.4 | 12 |
| Laguna del Hunco | 52 | 96.2 | 97.4 | 71 |
| Republic | 51.18 | 81.6 | 86.6 | 17 |
| Bonanza | 47.3 | 113.6 | 117.6 | 16 |

The PIP cross-validation RMSE was 0.136 on the log₁₀(LMA) scale across 107 extant sites.

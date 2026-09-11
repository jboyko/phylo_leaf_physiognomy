# PART 4: Empirical Application

## Methods

### Fossil leaf data

Leaf morphological measurements were taken from 10 late Cretaceous through early Eocene fossil sites from Peppe et al. (2011), spanning ages from 66.5 Ma (Fox Hills) to 47.3 Ma (Bonanza). Palacio de los Loros comprises the PL1 and PL2 collections, combined before species means are calculated. The dataset contains 360 species-by-site rows, with one row per fossil species at each site at which it occurs. Each row records the same 12 fossil-measurable predictors used to train the climate models described in Part 1. Each fossil species entry also carries four taxonomic fields (species, genus, family, order) and the site age in millions of years (age_ma). These fields are required for phylogenetic placement and are the same columns expected by the `dilp_pgls()` function in the dilp R package. A subset of fossil species also carry measurements of PW²/A (petiole width squared divided by blade area), the single predictor used in the LMA model; species without this measurement were excluded from LMA prediction.

### Phylogenetic placement

Each fossil occurrence was grafted onto the angiosperm scaffold tree (described in Part 1) using a genus-then-family-then-order fallback hierarchy and that occurrence's site age. The formal-only analysis is primary: a genus, family, or order reported in quotation marks is treated as informal and is not used as placement evidence. A separate sensitivity scenario provisionally includes those quoted ranks. If two or more extant scaffold tips match the highest available formal rank, their MRCA is the placement target; a single matching tip defines the target branch; otherwise the next formal rank is used, with the root as the final fallback. The edge length is chosen so that the fossil tip falls at its correct time depth. When a fossil predates the crown age of its placement clade, the code walks up the tree to find a branch alive at the fossil's age and splits that branch. This ancestral-branch approach ensures that no fossil is placed at a node younger than itself.

### PIP prediction

After placement, the phylogenetic variance-covariance matrix was computed for the combined set of extant training species and fossil occurrences using `vcv()` applied to the pruned tree. The cross-covariance block $V_{cross}$ between the $n$ training species and the $m$ fossil occurrences was extracted and scaled by Pagel's $\lambda$ estimated during PGLS fitting. The phylogenetic adjustment $V_{cross}^T V_{inv} e$ was computed as described in Part 1, with $V_{inv}$ and $e$ taken from the fitted PGLS components stored during model training. The trait-based component $X\beta$ was computed using site-specific trait means for each fossil occurrence. Each occurrence is also placed at its own age, so both its cross-covariance and phylogenetic adjustment can differ when the same label occurs at more than one site or age. Missing trait values were filled by bagged-tree imputation before constructing the design matrix. Predictions on the log(MAP) scale were averaged within each site and then exponentiated, giving geometric site means. Site-level MAT estimates remain arithmetic means of the occurrence-level predictions.

LMA predictions followed the same PIP framework using a separate model fitted to the extant LMA calibration dataset. The single predictor is log₁₀(PW²/A) (Dana Royer pers. comm.). For each fossil species with a valid PW²/A measurement, the design matrix $X$ contains an intercept column and a column of log₁₀(PW²/A) values aggregated to within-site species means. No imputation was applied; species without PW²/A measurements were excluded. Each LMA occurrence is placed at its own site age, with $V_{cross}$, $V_{inv}$, and $e$ drawn from the LMA model’s fitted components. Predictions on the log₁₀(LMA) scale were back-transformed as $10^{\hat{y}}$ before site averaging.

The climate predictions can be reproduced using the updated local `dilp_pgls()` implementation; the previously pinned public commit does not contain these changes. The function accepts a specimen-level data frame with the standard DiLP trait columns plus the species, genus, family, order, and age_ma fields, and returns site-level and species-level predictions together with a placement log recording how each fossil species was grafted onto the scaffold phylogeny.

## Results

All 360 fossil occurrences were retained in the formal-only climate analysis. PIP MAT estimates range from 14.0 °C at Republic to 21.6 °C at Cerrejon, and MAP estimates from 146 to 217 cm at the same sites. The practical non-phylogenetic comparison uses 12-trait site regression with imputation. Both methods use imputation for MAT and MAP.

**Table 6.** Site-level climate estimates for 10 fossil floras. PIP uses occurrence-specific traits and ages with formal-only taxonomy; LM is the 12-trait site regression with imputation. n is the number of species-site occurrences contributing to each PIP estimate. PIP MAP values are geometric means in cm; LM predicts directly at site level.

| Site | Age (Ma) | MAT PIP (°C) | MAT LM (°C) | MAP PIP (cm) | MAP LM (cm) | n |
| --- | --- | --- | --- | --- | --- | --- |
| Fox Hills | 66.5 | 16.2 | 18.1 | 152 | 83 | 24 |
| Williston Basin I | 64.75 | 17.5 | 14.6 | 168 | 137 | 20 |
| Williston Basin II | 63.5 | 16.6 | 13.8 | 161 | 101 | 23 |
| Palacio de los Loros | 61.7 | 17.2 | 13.5 | 166 | 107 | 29 |
| Williston Basin III | 59.75 | 16.6 | 16.1 | 160 | 106 | 18 |
| Cerrejon | 58 | 21.6 | 25.9 | 217 | 286 | 45 |
| Hubble Bubble | 55.8 | 19.2 | 17.8 | 167 | 96 | 16 |
| Laguna del Hunco | 51.9 | 16.6 | 11 | 168 | 139 | 119 |
| Republic | 49.4 | 14 | 8.7 | 146 | 66 | 41 |
| Bonanza | 47.3 | 17.1 | 10.3 | 148 | 116 | 25 |

Cross-validation RMSE is 3.414 °C for imputed PIP MAT and 0.525 for ln(MAP), the latter corresponding to a multiplicative error scale of approximately ×1.69. These are calibration error benchmarks, not fossil-specific confidence intervals, and do not include taxonomic-placement or age uncertainty. Under primary formal-only placement, 192 of the 360 occurrences attach at the root, 99 at family, 49 at genus, and 20 at order level. Provisionally including quoted taxonomic ranks changes rounded site PIP estimates by at most 0.5 °C and 3 cm MAP, both at Cerrejon. The sensitivity results do not resolve the validity of the provisional taxonomic assignments.


PIP LMA estimates range from 74.6 g m⁻² at Williston Basin III to 125.5 g m⁻² at Bonanza. The analysis includes 192 species-by-site occurrences across 9 sites with usable PW²/A measurements. Palacio de los Loros contributes 23 occurrences. Site-level estimates are given in Table 7.

**Table 7.** PIP and LM site-level LMA estimates for 9 fossil sites from Peppe et al. (2011). LMA is in g m⁻². n is the number of fossil species with PW²/A measurements contributing to the site mean.

| Site | Age (Ma) | LMA PIP (g m⁻²) | LMA LM (g m⁻²) | n |
| --- | --- | --- | --- | --- |
| Fox Hills | 66.5 | 89.1 | 86.9 | 16 |
| Williston Basin I | 64.75 | 95.9 | 92.4 | 11 |
| Williston Basin II | 63.5 | 80.2 | 79.2 | 13 |
| Palacio de los Loros | 61.7 | 86 | 83.9 | 23 |
| Williston Basin III | 59.75 | 74.6 | 73.2 | 13 |
| Hubble Bubble | 55.8 | 93 | 99.3 | 12 |
| Laguna del Hunco | 51.9 | 98.7 | 100.3 | 71 |
| Republic | 49.4 | 82.7 | 87.1 | 17 |
| Bonanza | 47.3 | 125.5 | 133.2 | 16 |

The cross-validation RMSE for LMA is 0.130 log₁₀ g m⁻², corresponding to a multiplicative uncertainty of approximately ×1.35 on the linear scale. Sites with few contributing species have additional uncertainty in their mean estimates.

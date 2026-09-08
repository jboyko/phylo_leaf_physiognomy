# PART 4: Empirical Application

## Methods

### Fossil leaf data

Leaf morphological measurements were taken from 11 late Cretaceous through early Eocene fossil sites from Peppe et al. (2011), spanning ages from 66.5 Ma (Fox Hills) to 47.3 Ma (Bonanza). The dataset contains 361 species-by-site rows, with one row per fossil species at each site at which it occurs. Each row records the same 12 fossil-measurable predictors used to train the climate models described in Part 1. Each fossil species entry also carries four taxonomic fields (species, genus, family, order) and the site age in millions of years (age_ma). These fields are required for phylogenetic placement and are the same columns expected by the `dilp_pgls()` function in the dilp R package. A subset of fossil species also carry measurements of PW²/A (petiole width squared divided by blade area), the single predictor used in the LMA model; species without this measurement were excluded from LMA prediction.

### Phylogenetic placement

Each fossil occurrence was grafted onto the angiosperm scaffold tree (described in Part 1) using a genus-then-family-then-order fallback hierarchy and that occurrence's site age. The formal-only analysis is primary: a genus, family, or order reported in quotation marks is treated as informal and is not used as placement evidence. A separate sensitivity scenario provisionally includes those quoted ranks. If two or more extant scaffold tips match the highest available formal rank, their MRCA is the placement target; a single matching tip defines the target branch; otherwise the next formal rank is used, with the root as the final fallback. The edge length is chosen so that the fossil tip falls at its correct time depth. When a fossil predates the crown age of its placement clade, the code walks up the tree to find a branch alive at the fossil's age and splits that branch. This ancestral-branch approach ensures that no fossil is placed at a node younger than itself.

### PIP prediction

After placement, the phylogenetic variance-covariance matrix was computed for the combined set of extant training species and fossil occurrences using `vcv()` applied to the pruned tree. The cross-covariance block $V_{cross}$ between the $n$ training species and the $m$ fossil occurrences was extracted and scaled by Pagel's $\lambda$ estimated during PGLS fitting. The phylogenetic adjustment $V_{cross}^T V_{inv} e$ was computed as described in Part 1, with $V_{inv}$ and $e$ taken from the fitted PGLS components stored during model training. The trait-based component $X\beta$ was computed using site-specific trait means for each fossil occurrence. Each occurrence is also placed at its own age, so both its cross-covariance and phylogenetic adjustment can differ when the same label occurs at more than one site or age. Missing trait values were filled by bagged-tree imputation before constructing the design matrix. Predictions on the log(MAP) scale were averaged within each site and then exponentiated, giving geometric site means. Site-level MAT estimates remain arithmetic means of the occurrence-level predictions.

LMA predictions followed the same PIP framework using a separate model fitted to the extant LMA calibration dataset. The single predictor is log₁₀(PW²/A) (Dana Royer pers. comm.). For each fossil species with a valid PW²/A measurement, the design matrix $X$ contains an intercept column and a column of log₁₀(PW²/A) values aggregated to within-site species means. No imputation was applied; species without PW²/A measurements were excluded. The separate LMA pipeline used species-level fossil placement at mean occurrence ages, with $V_{cross}$, $V_{inv}$, and $e$ drawn from its own fitted components. The occurrence-specific climate placement update described above has not been applied to these LMA results. Predictions on the log₁₀(LMA) scale were back-transformed as $10^{\hat{y}}$ before site averaging.

The climate predictions can be reproduced using the updated local `dilp_pgls()` implementation; the previously pinned public commit does not contain these changes. The function accepts a specimen-level data frame with the standard DiLP trait columns plus the species, genus, family, order, and age_ma fields, and returns site-level and species-level predictions together with a placement log recording how each fossil species was grafted onto the scaffold phylogeny.

## Results

All 361 fossil occurrences were retained in the formal-only climate analysis. PIP MAT estimates range from 14.1 °C at Republic to 21.6 °C at Cerrejon, and MAP estimates from 148 to 220 cm at the same sites. The practical non-phylogenetic comparison uses LM site sp+zero with imputation. Both methods use imputation for MAT and MAP.

**Table 6.** Site-level climate estimates for 11 fossil floras. PIP uses occurrence-specific traits and ages with formal-only taxonomy; LM is the site sp+zero imputed baseline. n is the number of species-site occurrences contributing to each PIP estimate. MAP values are geometric means in cm.

| Site | Age (Ma) | MAT PIP (°C) | MAT LM (°C) | MAP PIP (cm) | MAP LM (cm) | n |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Fox Hills | 66.5 | 16.3 | 18.1 | 153 | 83 | 24 |
| Williston Basin I | 64.75 | 17.5 | 14.6 | 168 | 138 | 20 |
| Williston Basin II | 63.5 | 16.7 | 13.8 | 161 | 102 | 23 |
| Palacio de los Loros PL1 | 61.7 | 17.1 | 12.9 | 168 | 118 | 24 |
| Palacio de los Loros PL2 | 61.7 | 18.2 | 17.6 | 161 | 83 | 6 |
| Williston Basin III | 59.75 | 16.6 | 16.1 | 160 | 106 | 18 |
| Cerrejon | 58 | 21.6 | 25.9 | 216 | 288 | 45 |
| Hubble Bubble | 55.8 | 19.2 | 17.8 | 166 | 96 | 16 |
| Laguna del Hunco | 51.9 | 16.7 | 11.0 | 168 | 139 | 119 |
| Republic | 49.4 | 14.1 | 8.6 | 146 | 67 | 41 |
| Bonanza | 47.3 | 17.1 | 10.3 | 149 | 116 | 25 |

Cross-validation RMSE is 3.414 °C for imputed PIP MAT and 0.525 for ln(MAP), the latter corresponding to a multiplicative error scale of approximately ×1.69. These are calibration error benchmarks, not fossil-specific confidence intervals, and do not include taxonomic-placement or age uncertainty. Under primary formal-only placement, 193 of the 361 occurrences attach at the root, 99 at family, 49 at genus, and 20 at order level. Provisionally including quoted taxonomic ranks changes rounded site PIP estimates by at most 0.5 °C and 4 cm MAP, both at Cerrejon. The sensitivity results do not resolve the validity of the provisional taxonomic assignments.


PIP LMA estimates ranged from 74.7 g m⁻² at Williston Basin III (59.75 Ma) to 124.6 g m⁻² at Bonanza (47.3 Ma). LM and PIP estimates were closely aligned across all sites, with differences of 10 g m⁻² or less except at Hubble Bubble and Bonanza. The number of species contributing to each site LMA average ranged from 4 at Palacio de los Loros PL2 to 71 at Laguna del Hunco. Site-level LMA predictions are given in Table 7.

**Table 7.** PIP and LM site-level LMA estimates for 11 fossil sites from Peppe et al. (2011). LMA is in g m⁻². n is the number of fossil species with PW²/A measurements contributing to the site mean.

| Site | Age (Ma) | LMA PIP (g m⁻²) | LMA LM (g m⁻²) | n |
| --- | --- | --- | --- | --- |
| Fox Hills | 66.5 | 89.1 | 86.9 | 153 |
| Williston Basin I | 64.8 | 96.8 | 93.0 | 168 |
| Williston Basin II | 63.5 | 80.5 | 79.2 | 161 |
| Palacio de los Loros PL1 | 61.7 | 79.3 | 79.7 | 168 |
| Palacio de los Loros PL2 | 61.7 | 95.2 | 100.3 | 161 |
| Williston Basin III | 59.8 | 74.7 | 73.2 | 160 |
| Cerrejon | 58.0 | 90.7 | 94.9 | 216 |
| Hubble Bubble | 55.8 | 92.8 | 99.3 | 166 |
| Laguna del Hunco | 51.9 | 98.1 | 100.3 | 168 |
| Republic | 49.4 | 83.4 | 87.2 | 146 |
| Bonanza | 47.3 | 124.6 | 133.2 | 149 |

The cross-validation RMSE for LMA is 0.130 log₁₀ g m⁻², corresponding to a multiplicative uncertainty of approximately ×1.35 on the linear scale. As with the climate estimates, sites with few contributing species carry additional uncertainty; Palacio de los Loros PL2 (n = 4) should be interpreted with particular caution.

# PART 4: Empirical Application

## Methods

### Fossil leaf data

Leaf morphological measurements were taken from 11 late Cretaceous through early Eocene fossil sites from Peppe et al. (2011), spanning ages from 66.5 Ma (Fox Hills) to 47.3 Ma (Bonanza). The dataset contains 361 species-by-site rows, with one row per fossil species at each site at which it occurs. Each row records the same 12 fossil-measurable predictors used to train the climate models described in Part 1. Each fossil species entry also carries four taxonomic fields (species, genus, family, order) and the site age in millions of years (age_ma). These fields are required for phylogenetic placement and are the same columns expected by the `dilp_pgls()` function in the dilp R package. A subset of fossil species also carry measurements of PW²/A (petiole width squared divided by blade area), the single predictor used in the LMA model; species without this measurement were excluded from LMA prediction.

### Phylogenetic placement

Each fossil species × site occurrence was grafted independently onto the angiosperm scaffold tree (described in Part 1) using that site's age and a genus-then-family-then-order fallback hierarchy. If one or more extant tips with a matching genus were present in the scaffold, the occurrence was attached at the most recent common ancestor (MRCA) of those tips. If no genus match was found, the MRCA of extant tips belonging to the same family was used, followed by order and then root. Placement targets were resolved only from the original extant scaffold tips; previously grafted fossils were never used as taxonomic evidence for later placements. The edge length assigned to each graft was chosen so that the fossil tip falls at the correct time depth, computed from the calibrated scaffold tree height minus the occurrence's site age in Ma. When a fossil predates the crown age of its placement clade, the code walks up the tree to find the branch that was alive at the fossil's age and splits that branch, attaching the fossil with a zero-length terminal edge. Thus, a species label occurring at multiple sites receives a separate, time-correct tip and phylogenetic adjustment at every site. The same rule is used when genus, family, or order is the finest available identification.

### PIP prediction

After placement, the phylogenetic variance-covariance matrix was computed for the combined set of extant training species and fossil occurrences using `vcv()` applied to the pruned tree. The cross-covariance block $V_{cross}$ between the $n$ training species and the $m$ fossil occurrences was extracted and scaled by Pagel's $\lambda$ estimated during PGLS fitting. The phylogenetic adjustment $V_{cross}^T V_{inv} e$ was computed as described in Part 1, with $V_{inv}$ and $e$ taken from the fitted PGLS components stored during model training. The trait-based component $X\beta$ was computed using species-within-site trait means. Consequently, species appearing at more than one site contribute distinct trait values, placement ages, and phylogenetic adjustments at each site; no fossil physiognomic trait or age is averaged across sites in the primary model. Missing trait values were filled by bagged-tree imputation before constructing the design matrix. Predictions on the log(MAP) scale were exponentiated before site averaging. Site-level MAT and MAP estimates are the mean of occurrence-level predictions within each site.

LMA predictions followed the same occurrence-specific PIP framework using a separate model fitted to the extant LMA calibration dataset. The single predictor is log₁₀(PW²/A) (Dana Royer pers. comm.). Each fossil species × site occurrence with a valid PW²/A measurement was placed at its own site age, and the design matrix $X$ contained an intercept plus its within-site species mean log₁₀(PW²/A). No imputation was applied; occurrences without PW²/A measurements were excluded. Predictions on the log₁₀(LMA) scale were back-transformed as $10^{\hat{y}}$ before site averaging.

The climate analysis can also be run through the `dilp_pgls()` function in the dilp R package. The function accepts a specimen-level data frame with the standard DiLP trait columns plus the species, genus, family, order, and age_ma fields, and returns site-level and species-by-site occurrence predictions together with a placement log recording how each fossil occurrence was grafted onto the scaffold phylogeny.

## Results

PIP MAT estimates ranged from 14.0 °C at Republic (49.4 Ma) to 21.6 °C at Cerrejon (58.0 Ma). MAP estimates ranged from 148 cm at Republic to 220 cm at Cerrejon. The number of fossil occurrences contributing to each site average ranged from 6 at Palacio de los Loros PL2 to 119 at Laguna del Hunco. Site-level predictions are given in Table 4.

**Table 4.** PIP site-level paleoclimate estimates for 11 fossil sites from Peppe et al. (2011). MAT is in °C. MAP is in cm. n is the number of fossil species × site occurrences contributing to the site mean.

| Site | Age (Ma) | MAT (°C) | MAP (cm) | n |
| --- | --- | --- | --- | --- |
| Fox Hills | 66.5 | 16.2 | 155 | 24 |
| Williston Basin I | 64.8 | 17.5 | 172 | 20 |
| Williston Basin II | 63.5 | 16.6 | 164 | 23 |
| Palacio de los Loros PL1 | 61.7 | 17.0 | 170 | 24 |
| Palacio de los Loros PL2 | 61.7 | 18.1 | 164 | 6 |
| Williston Basin III | 59.8 | 16.6 | 162 | 18 |
| Cerrejon | 58.0 | 21.6 | 220 | 45 |
| Hubble Bubble | 55.8 | 19.2 | 170 | 16 |
| Laguna del Hunco | 51.9 | 16.6 | 169 | 119 |
| Republic | 49.4 | 14.0 | 148 | 41 |
| Bonanza | 47.3 | 17.1 | 150 | 25 |

The cross-validation RMSE from Part 3 provides the relevant uncertainty benchmark for these estimates. PIP MAT predictions carry an expected error of approximately 3.4 °C (RMSE across 93 held-out modern sites). log(MAP) RMSE of 0.52 log cm corresponds to a multiplicative uncertainty of roughly a factor of 1.7 on the linear precipitation scale, so MAP estimates should be treated as order-of-magnitude reconstructions. Sites with few contributing occurrences, particularly Palacio de los Loros PL2 (n = 6), carry additional uncertainty because the site mean is based on a small sample of occurrence-level predictions.

PIP LMA estimates ranged from 74.6 g m⁻² at Williston Basin III (59.75 Ma) to 125.5 g m⁻² at Bonanza (47.3 Ma). The number of occurrences contributing to each site LMA average ranged from 4 at Palacio de los Loros PL2 to 71 at Laguna del Hunco. Cerrejon was omitted because the current fossil dataset contains no PW²/A measurements for that site. Site-level LMA predictions are given in Table 5.

**Table 5.** PIP and LM site-level LMA estimates for the 10 fossil sites with PW²/A data. LMA is in g m⁻². n is the number of fossil species × site occurrences with PW²/A measurements contributing to the site mean.

| Site | Age (Ma) | LMA PIP (g m⁻²) | LMA LM (g m⁻²) | n |
| --- | --- | --- | --- | --- |
| Fox Hills | 66.5 | 89.1 | 86.9 | 16 |
| Williston Basin I | 64.8 | 95.9 | 92.4 | 11 |
| Williston Basin II | 63.5 | 80.2 | 79.2 | 13 |
| Palacio de los Loros PL1 | 61.7 | 82.7 | 80.5 | 20 |
| Palacio de los Loros PL2 | 61.7 | 100.5 | 100.3 | 4 |
| Williston Basin III | 59.8 | 74.6 | 73.2 | 13 |
| Hubble Bubble | 55.8 | 93.0 | 99.3 | 12 |
| Laguna del Hunco | 51.9 | 98.7 | 100.3 | 71 |
| Republic | 49.4 | 82.7 | 87.1 | 17 |
| Bonanza | 47.3 | 125.5 | 133.2 | 16 |

The cross-validation RMSE for LMA is 0.130 log₁₀ g m⁻², corresponding to a multiplicative uncertainty of approximately ×1.35 on the linear scale. As with the climate estimates, sites with few contributing species carry additional uncertainty; Palacio de los Loros PL2 (n = 4) should be interpreted with particular caution.

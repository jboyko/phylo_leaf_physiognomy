# PART 4: Empirical Application

## Methods

### Fossil leaf data

Leaf morphological measurements were taken from 11 late Cretaceous through early Eocene fossil sites from Peppe et al. (2011), spanning ages from 66.5 Ma (Fox Hills) to 47.3 Ma (Bonanza). The dataset contains 361 species-by-site rows, with one row per fossil species at each site at which it occurs. Each row records the same 12 fossil-measurable predictors used to train the climate models described in Part 1. Each fossil species entry also carries four taxonomic fields (species, genus, family, order) and the site age in millions of years (age_ma). These fields are required for phylogenetic placement and are the same columns expected by the `dilp_pgls()` function in the dilp R package. A subset of fossil species also carry measurements of PW²/A (petiole width squared divided by blade area), the single predictor used in the LMA model; species without this measurement were excluded from LMA prediction.

### Phylogenetic placement

Each fossil species was grafted onto the angiosperm scaffold tree (described in Part 1) using a genus-then-family-then-order fallback hierarchy. Placement proceeded as follows. If one or more extant tips with a matching genus were present in the scaffold, the fossil was attached at the most recent common ancestor (MRCA) of those tips. If no genus match was found, the MRCA of tips belonging to the same family was used. If no family match was found, the MRCA of tips belonging to the same order was used. If no order match was found, the fossil was attached at the root with a warning. The edge length assigned to each graft was chosen so that the fossil tip falls at the correct time depth, computed from the calibrated scaffold tree height minus the fossil's age in Ma. When a fossil predates the crown age of its placement clade, the code walks up the tree to find the branch that was alive at the fossil's age and splits that branch, attaching the fossil with a zero-length terminal edge. This ancestral branch approach ensures that no fossil is placed at a node younger than itself.

### PIP prediction

After placement, the phylogenetic variance-covariance matrix was computed for the combined set of extant training species and fossil species using `vcv()` applied to the pruned tree. The cross-covariance block $V_{cross}$ between the $n$ training species and the $m$ fossil species was extracted and scaled by Pagel's $\lambda$ estimated during PGLS fitting. The phylogenetic adjustment $V_{cross}^T V_{inv} e$ was computed as described in Part 1, with $V_{inv}$ and $e$ taken from the fitted PGLS components stored during model training. The trait-based component $X\beta$ was computed using site-specific trait means for each fossil species, so that species appearing at more than one site contribute a distinct $X\beta$ value for each site while sharing the same phylogenetic adjustment. Missing trait values were filled by bagged-tree imputation before constructing the design matrix. Predictions on the log(MAP) scale were exponentiated before site averaging. Site-level MAT and MAP estimates are the mean of species-level predictions within each site.

LMA predictions followed the same PIP framework using a separate model fitted to the extant LMA calibration dataset. The single predictor is log₁₀(PW²/A) (Dana Royer pers. comm.). For each fossil species with a valid PW²/A measurement, the design matrix $X$ contains an intercept column and a column of log₁₀(PW²/A) values aggregated to within-site species means. No imputation was applied; species without PW²/A measurements were excluded. The phylogenetic adjustment used the same grafted tree and cross-covariance procedure described above, with $V_{cross}$, $V_{inv}$, and $e$ drawn from the LMA model components. Predictions on the log₁₀(LMA) scale were back-transformed as $10^{\hat{y}}$ before site averaging.

These analyses can be reproduced using the `dilp_pgls()` function in the dilp R package. The function accepts a specimen-level data frame with the standard DiLP trait columns plus the species, genus, family, order, and age_ma fields, and returns site-level and species-level predictions together with a placement log recording how each fossil species was grafted onto the scaffold phylogeny.

## Results

PIP MAT estimates ranged from 14.2 °C at Republic (49.4 Ma) to 22.1 °C at Cerrejon (58.0 Ma). MAP estimates ranged from 151 cm at Republic and Bonanza to 225 cm at Cerrejon. The number of species contributing to each site average ranged from 6 at Palacio de los Loros PL2 to 118 at Laguna del Hunco. Site-level predictions are given in Table 4.

**Table 4.** PIP site-level paleoclimate estimates for 11 fossil sites from Peppe et al. (2011). MAT is in °C. MAP is in cm. n is the number of fossil species contributing to the site mean.

| Site | Age (Ma) | MAT (°C) | MAP (cm) | n |
| --- | --- | --- | --- | --- |
| Fox Hills | 66.5 | 16.2 | 156 | 25 |
| Williston Basin I | 64.8 | 17.1 | 172 | 20 |
| Williston Basin II | 63.5 | 16.2 | 166 | 23 |
| Palacio de los Loros PL1 | 61.7 | 17.0 | 165 | 26 |
| Palacio de los Loros PL2 | 61.7 | 18.1 | 159 | 6 |
| Williston Basin III | 59.8 | 16.3 | 163 | 18 |
| Cerrejon | 58.0 | 22.1 | 225 | 45 |
| Hubble Bubble | 55.8 | 19.8 | 170 | 16 |
| Laguna del Hunco | 51.9 | 16.8 | 170 | 118 |
| Republic | 49.4 | 14.2 | 151 | 39 |
| Bonanza | 47.3 | 17.0 | 151 | 25 |

The cross-validation RMSE from Part 3 provides the relevant uncertainty benchmark for these estimates. PIP MAT predictions carry an expected error of approximately 3.4 °C (RMSE across 93 held-out modern sites). log(MAP) RMSE of 0.52 log cm corresponds to a multiplicative uncertainty of roughly a factor of 1.7 on the linear precipitation scale, so MAP estimates should be treated as order-of-magnitude reconstructions. Sites with few contributing species, particularly Palacio de los Loros PL2 (n = 6), carry additional uncertainty because the site mean is based on a small sample of species-level predictions.

PIP LMA estimates ranged from 74.7 g m⁻² at Williston Basin III (59.75 Ma) to 124.6 g m⁻² at Bonanza (47.3 Ma). LM and PIP estimates were closely aligned across all sites, with differences of 10 g m⁻² or less except at Hubble Bubble and Bonanza. The number of species contributing to each site LMA average ranged from 4 at Palacio de los Loros PL2 to 71 at Laguna del Hunco. Site-level LMA predictions are given in Table 5.

**Table 5.** PIP and LM site-level LMA estimates for 11 fossil sites from Peppe et al. (2011). LMA is in g m⁻². n is the number of fossil species with PW²/A measurements contributing to the site mean.

| Site | Age (Ma) | LMA PIP (g m⁻²) | LMA LM (g m⁻²) | n |
| --- | --- | --- | --- | --- |
| Fox Hills | 66.5 | 89.1 | 86.9 | 16 |
| Williston Basin I | 64.8 | 96.8 | 93.0 | 11 |
| Williston Basin II | 63.5 | 80.5 | 79.2 | 13 |
| Palacio de los Loros PL1 | 61.7 | 79.3 | 79.7 | 22 |
| Palacio de los Loros PL2 | 61.7 | 95.2 | 100.3 | 4 |
| Williston Basin III | 59.8 | 74.7 | 73.2 | 13 |
| Cerrejon | 58.0 | 90.7 | 94.9 | 24 |
| Hubble Bubble | 55.8 | 92.8 | 99.3 | 12 |
| Laguna del Hunco | 51.9 | 98.1 | 100.3 | 71 |
| Republic | 49.4 | 83.4 | 87.2 | 17 |
| Bonanza | 47.3 | 124.6 | 133.2 | 16 |

The cross-validation RMSE for LMA is 0.130 log₁₀ g m⁻², corresponding to a multiplicative uncertainty of approximately ×1.35 on the linear scale. As with the climate estimates, sites with few contributing species carry additional uncertainty; Palacio de los Loros PL2 (n = 4) should be interpreted with particular caution.

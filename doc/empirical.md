# PART 4: Empirical Application

## Methods

### Fossil leaf data

Leaf morphological measurements were taken from 11 late Cretaceous through early Eocene fossil sites from Peppe et al. (2011), spanning ages from 66.5 Ma (Fox Hills) to 47.3 Ma (Bonanza). The dataset contains 361 species-by-site rows, with one row per fossil species at each site at which it occurs. Each row records the same 12 fossil-measurable predictors used to train the models described in Part 1. Each fossil species entry also carries four taxonomic fields (species, genus, family, order) and the site age in millions of years (age_ma). These fields are required for phylogenetic placement and are the same columns expected by the `dilp_pgls()` function in the dilp R package.

### Phylogenetic placement

Each fossil species was grafted onto the angiosperm scaffold tree (described in Part 1) using a genus-then-family-then-order fallback hierarchy. Placement proceeded as follows. If one or more extant tips with a matching genus were present in the scaffold, the fossil was attached at the most recent common ancestor (MRCA) of those tips. If no genus match was found, the MRCA of tips belonging to the same family was used. If no family match was found, the MRCA of tips belonging to the same order was used. If no order match was found, the fossil was attached at the root with a warning. The edge length assigned to each graft was chosen so that the fossil tip falls at the correct time depth, computed from the calibrated scaffold tree height minus the fossil's age in Ma. When a fossil predates the crown age of its placement clade, the code walks up the tree to find the branch that was alive at the fossil's age and splits that branch, attaching the fossil with a zero-length terminal edge. This ancestral branch approach ensures that no fossil is placed at a node younger than itself.

### PIP prediction

After placement, the phylogenetic variance-covariance matrix was computed for the combined set of extant training species and fossil species using `vcv()` applied to the pruned tree. The cross-covariance block $V_{cross}$ between the $n$ training species and the $m$ fossil species was extracted and scaled by Pagel's $\lambda$ estimated during PGLS fitting. The phylogenetic adjustment $V_{cross}^T V_{inv} e$ was computed as described in Part 1, with $V_{inv}$ and $e$ taken from the fitted PGLS components stored during model training. The trait-based component $X\beta$ was computed using site-specific trait means for each fossil species, so that species appearing at more than one site contribute a distinct $X\beta$ value for each site while sharing the same phylogenetic adjustment. Missing trait values were filled by bagged-tree imputation before constructing the design matrix. Predictions on the log(MAP) scale were exponentiated before site averaging. Site-level MAT and MAP estimates are the mean of species-level predictions within each site.

These analyses can be reproduced using the `dilp_pgls()` function in the dilp R package. The function accepts a specimen-level data frame with the standard DiLP trait columns plus the species, genus, family, order, and age_ma fields, and returns site-level and species-level predictions together with a placement log recording how each fossil species was grafted onto the scaffold phylogeny.

## Results

PIP MAT estimates ranged from 14.4 °C at Republic (49.4 Ma) to 22.6 °C at Cerrejon (58.0 Ma). MAP estimates ranged from 149 cm at Bonanza (47.3 Ma) to 222 cm at Cerrejon. The number of species contributing to each site average ranged from 6 at Palacio de los Loros PL2 to 118 at Laguna del Hunco. Site-level predictions are given in Table 4.

**Table 4.** PIP site-level paleoclimate estimates for 11 fossil sites from Peppe et al. (2011). MAT is in °C. MAP is in cm. n is the number of fossil species contributing to the site mean.

| Site | Age (Ma) | MAT (°C) | MAP (cm) | n |
| --- | --- | --- | --- | --- |
| Fox Hills | 66.5 | 16.5 | 157 | 25 |
| Williston Basin I | 64.8 | 17.6 | 172 | 20 |
| Williston Basin II | 63.5 | 16.6 | 170 | 23 |
| Palacio de los Loros PL1 | 61.7 | 17.5 | 171 | 26 |
| Palacio de los Loros PL2 | 61.7 | 18.5 | 162 | 6 |
| Williston Basin III | 59.8 | 16.8 | 166 | 18 |
| Cerrejon | 58.0 | 22.6 | 222 | 45 |
| Hubble Bubble | 55.8 | 20.4 | 168 | 16 |
| Laguna del Hunco | 51.9 | 17.3 | 173 | 118 |
| Republic | 49.4 | 14.4 | 153 | 39 |
| Bonanza | 47.3 | 17.7 | 149 | 25 |

The cross-validation RMSE from Part 3 provides the relevant uncertainty benchmark for these estimates. PIP MAT predictions carry an expected error of approximately 3.6 °C (RMSE across 92 held-out modern sites). log(MAP) RMSE of 0.54 log cm corresponds to a multiplicative uncertainty of roughly a factor of 1.7 on the linear precipitation scale, so MAP estimates should be treated as order-of-magnitude reconstructions. Sites with few contributing species, particularly Palacio de los Loros PL2 (n = 6), carry additional uncertainty because the site mean is based on a small sample of species-level predictions.

# Climate model update for Dana

The updated MAT and MAP analysis is ready for review ahead of the October talk. PIP has lower cross-validated error than the requested non-phylogenetic site LM baseline. Fossil placement now preserves each occurrence’s site age, and all 361 occurrences contribute to the primary analysis.

## Calibration comparison

| Model | MAT RMSE (°C) | ln(MAP) RMSE | n sites MAT / MAP |
| --- | ---: | ---: | ---: |
| PIP impute | 3.414 | 0.525 | 93 |
| LM site sp+zero impute | 4.217 | 0.672 | 93 |
| Published DiLP in-sample | 3.791 | 0.554 | 93 / 91 |

PIP and LM use 10-fold cross-validation with whole sites held out and models and imputers refitted within each fold. Published DiLP uses fixed published coefficients scored on the calibration data; it is an in-sample historical reference. Its MAP score must not be presented as a like-for-like validation result. PIP versus the site LM is a practical comparison; PGLS versus PIP isolates the prediction-time phylogenetic correction.

MAP validation and fossil prediction use geometric site means: average species log predictions, then exponentiate for MAP in cm. On the same held-out sites, geometric averaging gave lower PIP impute RMSE than arithmetic averaging on both the log scale (0.525 versus 0.559) and the original scale (76.9 versus 81.2 cm). MAT is unchanged. The complete-case-training MAP variant has RMSE 0.508; this talk comparison retains the requested imputed variant.

## Fossil climate comparison

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

PIP uses formal-only taxonomy and each occurrence’s own site traits and age. Quoted ranks are excluded as placement evidence; including them provisionally changes site estimates by at most 0.5 °C and 4 cm MAP. Taxonomic anchors are restricted to the original extant scaffold, and all 361 occurrences are retained. The table corresponds to Table 6 after the validation section’s Tables 1–5.

## Interpretation for the talk

Fitting still averages extant traits across sites within each operational taxon, losing within-species variation. Prediction preserves site-specific traits and ages. Calibration labels with no species epithet can pool unnamed morphotypes and remain subject to taxonomic review. The benchmark tests new-site prediction with known extant relationships; it does not quantify fossil-placement uncertainty. Calibration RMSE is an error benchmark, not a fossil-specific confidence interval.

The requested comparison concerns MAT and MAP. LMA remains in the manuscript and has not been rerun for this update. Brian’s involvement is settled.

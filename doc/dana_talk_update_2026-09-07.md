# Climate model update for Dana

We compared three methods for estimating mean annual temperature (MAT) and mean annual precipitation (MAP) from leaf morphology. Phylogenetically informed prediction (PIP) had the lowest cross-validated error for both variables. We applied PIP and a 12-trait site regression to 10 fossil floras comprising 360 species-by-site occurrences.

## Models and validation

PIP is fitted to extant species’ mean leaf traits and climate. It predicts climate for each species at a new site using that species’ local traits, then adds a correction based on its relatedness to extant calibration species and their regression residuals. The species predictions are averaged to estimate the site’s climate.

The 12-trait site regression is fitted to site averages of 12 traits measurable on fossil leaves. Measurements are first averaged within species at each site, then across species to obtain site means. The model predicts climate directly from those site means.

For both methods, a leaf identified as untoothed receives zero for missing tooth-count and tooth-area traits, and one for a missing perimeter ratio, before traits are averaged. These values represent the absence of teeth. Other missing traits are estimated from observed traits using bagged regression trees fitted to the calibration data.

DiLP regression uses the three predictors and trait processing specified by Peppe et al. (2011) for each climate variable. MAT predictors are the percentage of untoothed species, feret diameter ratio, and tooth count per internal perimeter. MAP predictors are the natural logarithms of leaf area, tooth count per internal perimeter, and perimeter ratio. Log transformations are applied before site averaging.

We evaluated all three methods using the same ten groups of sites. In each round, one group was held out, coefficients were fitted using the remaining sites, and climate was predicted for the held-out sites. Imputation models were also fitted using only training data. Each site received one held-out prediction. RMSE summarizes the differences between predicted and observed climate; smaller values indicate better prediction.

## Calibration comparison

| Model | MAT RMSE (°C) | ln(MAP) RMSE | n sites MAT / MAP |
| --- | --- | --- | --- |
| PIP | 3.451 | 0.535 | 92 / 90 |
| 12-trait site regression | 4.269 | 0.695 | 92 / 90 |
| DiLP regression | 3.891 | 0.575 | 92 / 90 |

The calibration contains 92 sites. Yasuni-ridgetop and Yasuni-upper slope are treated as one site, while Yasuni-valley bottom remains separate. The MAT comparison includes all 92 sites; the MAP comparison includes the same 90 sites for every model. Kepong and Tanjung Tuan lack the logged tooth measurements needed for DiLP’s MAP equation. MAP errors are measured on the natural-log scale, with precipitation expressed in centimetres.

## Fossil climate comparison

Palacio de los Loros combines the PL1 and PL2 collections at an age of 64.08 Ma. Leaf measurements are averaged within each species at each fossil site. Each occurrence is placed on the extant phylogeny at its site’s age using its formal genus, family, or order identification. PIP combines that occurrence’s local traits with its phylogenetic correction. The 12-trait site regression predicts directly from the site’s mean traits. Missing fossil traits are estimated using the extant calibration imputers.

In the table, LM denotes the 12-trait site regression, age is in millions of years, and n is the number of species-by-site occurrences contributing to PIP. Site MAT estimates are arithmetic means of species predictions. Site MAP estimates are geometric means, calculated by averaging species predictions on the log scale and exponentiating.

| Site | Age (Ma) | MAT PIP (°C) | MAT LM (°C) | MAP PIP (cm) | MAP LM (cm) | n |
| --- | --- | --- | --- | --- | --- | --- |
| Fox Hills | 67 | 16.3 | 18.1 | 153 | 83 | 24 |
| Williston Basin I | 65.07 | 17.5 | 14.5 | 168 | 137 | 20 |
| Palacio de los Loros | 64.08 | 17.2 | 13.5 | 166 | 108 | 29 |
| Williston Basin II | 63.8 | 16.7 | 13.7 | 161 | 102 | 23 |
| Williston Basin III | 60.4 | 16.6 | 16 | 160 | 107 | 18 |
| Cerrejon | 59 | 21.6 | 25.9 | 216 | 283 | 45 |
| Hubble Bubble | 55.8 | 19.2 | 17.8 | 166 | 96 | 16 |
| Laguna del Hunco | 52 | 16.7 | 11 | 168 | 139 | 119 |
| Republic | 51.18 | 14.1 | 8.6 | 146 | 67 | 41 |
| Bonanza | 47.3 | 17.1 | 10.3 | 149 | 116 | 25 |

A taxonomic sensitivity analysis also uses provisional identifications reported in quotation marks. These placements change site PIP estimates by at most 0.5 °C for MAT and 4 cm for MAP.

## Interpretation

PIP predicts held-out extant sites more accurately than either site regression in this comparison. Its fossil estimates range from 14.1 to 21.6 °C for MAT and 146 to 216 cm for MAP.

Species-level calibration averages traits and climate across extant occurrences, so within-species variation is reduced during fitting. Some genus-level calibration labels combine unnamed morphotypes and require taxonomic review. Applying the validation results to fossils also depends on the accuracy of taxonomic assignments and fossil ages.

## Reproducibility

The site folds and held-out PIP and 12-trait regression predictions are stored in `tables/loso_cv_site_predictions.csv`. Running `Rscript code/03c_dilp_cv.R` fits DiLP on those folds and writes its predictions and coefficients to `tables/dilp_cv_site_predictions.csv` and `tables/dilp_cv_coefficients.csv`. The matching-site RMSE comparison is saved in `tables/dana_cv_comparison.csv`.

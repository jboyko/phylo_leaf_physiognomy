# Data S1

The `data/` folder contains the inputs for the PhyloDiLP analysis.

- `Peppe_2011_calibration_data_leaf_level_clean.csv` contains extant leaf measurements, updated taxonomy, and climate. MAT is in °C, source MAP is in mm, and LMA is in g/m². The analysis converts MAP to cm. [Peppe et al. (2011)](https://doi.org/10.1111/j.1469-8137.2010.03615.x), [Lowe et al. (2025)](https://doi.org/10.1002/ajb2.70019)
- `Peppe_2011_fossil_data_June_2026_leaf_level_clean.csv` contains Dana Royer’s June 2026 measurements and taxonomy for the Peppe fossil specimens, with the LM17 edits below.
- `extra_calibration_data_for_LMA.csv` adds extant species-by-site LMA records in g/m². It preserves the checked-in taxonomy. [Royer et al. (2007)](https://doi.org/10.1666/07001.1)
- `best_wcvp.tre_dated` is the project’s 123,182-tip WCVP tree snapshot, with branch lengths in millions of years. [Carruthers et al. (2026)](https://doi.org/10.1002/aps3.70092)
- `fossil_site_ages.csv` lists the analytical ages in Ma supplied by Dana Royer and Dan Peppe, as recorded in `code/fossil_taxonomy.R`. Palacio de los Loros PL1 and PL2 form one site at 64.08 Ma.
- `published_dilp_fossil.csv` gives the published DiLP estimates in °C and cm. Bonanza is `NA` for both. [Peppe et al. (2011), Table 1](https://doi.org/10.1111/j.1469-8137.2010.03615.x)

Leaf lengths are in cm, areas in cm², and tooth measurements are counts.

At Palacio de los Loros, the combined *Sterculia acuminataloba* occurrence uses Malvaceae and Malvales. For LM17 (*Averhoites affinus*), `Oxidales` was corrected to `Oxalidales`. Both its order and family are tentative. Quotation marks and trailing question marks flag ranks for use only in the sensitivity analysis. [Issue #39](https://github.com/jboyko/phylo_leaf_physiognomy/issues/39)

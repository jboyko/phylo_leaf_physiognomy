# Dana project update — September 8, 2026

The climate analysis repairs, refreshed validation and fossil results, local package update, and revised manuscript extracts are complete under the existing operational-taxon definitions. Brian's contact is settled. The main scientific decision still needing James and Dana's review is how unnamed calibration taxa should be grouped.

## Completed work

- Fossil climate predictions use each occurrence's site age and within-site traits. Taxonomic placement uses an immutable extant scaffold, preventing previously inserted fossils from becoming taxonomic anchors. Palacio de los Loros PL1 and PL2 are combined before species-site averaging. Both previously failed records now place successfully: Cornus nebrascensis at Williston Basin I and Crataegus sp. (rp42) at Republic. All 360 combined occurrences are retained.
- Following the averaging comparison requested by James, precipitation validation and fossil predictions now use geometric site means. The complete ten-fold climate CV was rerun and figures and diagnostics refreshed. MAT predictions were unchanged.
- The climate analysis and release code were synchronized. The local `/Users/jboyko/dilp` source now implements occurrence-specific predictions, target-specific trait imputation, formal versus informal taxonomy scenarios, explicit input validation, and visible placement failures. Exported package model data include the refreshed validation errors.
- Table 1, Table 6, methods, discussion, and a short talk handoff were updated in Markdown and Word. All pages of the five final Word documents were rendered and visually inspected. LMA results remain in the manuscript; its legacy prediction workflow was not rerun or modernized during this climate update.

## Current calibration comparison

| Configuration | MAT RMSE (°C) | ln(MAP) RMSE | Evaluation |
| --- | ---: | ---: | --- |
| PIP impute | 3.414 | 0.525 | Ten-fold site-grouped CV |
| LM site sp+zero impute | 4.217 | 0.672 | Same CV |
| Published DiLP coefficients | 3.791 | 0.554 | In-sample reference |

Geometric averaging was selected after comparison on the same 93 held-out sites: log RMSE 0.525 versus 0.559 for arithmetic averaging, original-scale RMSE 76.9 versus 81.2 cm, and log MAE 0.409 versus 0.424. The comparison is saved in `tables/map_averaging_comparison.csv`. This was a choice based on the existing CV results, not an independent validation of the selected rule. Archived geometric CV predictions were restored without model refitting; the two runs had identical folds, observations, and MAT predictions. Fossil site MAP values fall by roughly 1–4 cm, with MAT unchanged. PIP complete-case training has MAP RMSE 0.508. The best imputed site LM is the untoothed-excluded configuration (3.722 °C and 0.583 ln(MAP)); the simple talk comparison is not a comparison exclusively against the best alternatives. PIP versus site LM also changes fitting and aggregation, so its performance difference cannot be attributed solely to phylogeny.

## Verification

- Site-aggregation regression checks cover arithmetic versus geometric MAP means, missing values, and numerical overflow.
- All 360 fossil occurrences were inserted in original, reversed, and shuffled order. Extant–fossil covariance agreed to numerical tolerance, and final tip ages matched their supplied ages. An age older than the scaffold root fails without modifying the tree.
- Both formal-only and informal-inclusive fossil scenarios completed. Their largest rounded site differences were 0.5 °C and 4 cm MAP, both at Cerrejon. Formal-only placement includes 193 root, 99 family, 49 genus, and 20 order placements; successful placement does not imply precise taxonomic information.
- Local package tests passed. Source and temporarily installed package predictions matched the main pipeline for all 360 occurrences in both taxonomy scenarios, including unrounded MAT and MAP values. Existing raw-data quality and deprecation warnings remain visible.
- An isolated release-code run reproduced the four checked fossil output files byte for byte. Main and release climate scripts were also compared as parsed R expressions. The repository whitespace check passed.

## Decision retained for review

There are 1,740 calibration training labels; 213 have blank species epithets and 79 of those occur at multiple sites. For example, `Salix_` pools 55 site–morphotype entries across 20 sites. The review inventory is `doc/calibration_unnamed_taxa_review.csv`. Existing identifiers were retained to avoid making a new taxonomic assumption. James and Dana should decide whether unnamed entries represent shared taxa or require separate site–morphotype identifiers. Changing this requires rebuilding the training data, refitting models, repeating CV and fossil predictions, and refreshing package data and tables.

Even for correctly identified shared taxa, fitting grand means across sites removes within-taxon climate variation. The revised methods distinguish this training choice from the within-site traits and occurrence ages preserved at prediction time. Fossil uncertainty bars are calibration-error benchmarks, not fossil-specific confidence intervals.

## Deliverables and release status

- `doc/dana_talk_update_2026-09-07.docx`: short MAT/MAP handoff with current tables and limitations.
- `doc/model_validation.docx`: calibration methods, Table 1, and diagnostics; retained LMA tables.
- `doc/empirical.docx`: fossil methods and updated climate Table 6; retained LMA Table 7.
- `doc/fossil_predict.docx`: prediction explanation and worked example.
- `doc/discussion.docx`: revised interpretation and limitations.

Matching Markdown files are the editable text sources. Changes are local and uncommitted. The updated `dilp` package has not been published; the old pinned public version is not the updated prediction implementation. Nothing was emailed or pushed. The earlier correspondence audit is a historical record of findings before these repairs, not the current status.

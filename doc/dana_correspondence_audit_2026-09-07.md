# Dana correspondence audit — September 7, 2026

Follow-up: James confirmed on September 7 that Brian has been contacted and is
set. That action is closed. The findings below record the pre-repair audit;
the [September 8 work report](dana_update_work_report_2026-09-08.md) records fixes and refreshed results.

## Scope and conclusion

Audited James’s pasted August 6–26 correspondence against the local analysis repository, stored outputs, manuscript Markdown, release copy, and local `/Users/jboyko/dilp` checkout. Two Terra agents independently reviewed manuscript and implementation. This is an inspection of code and existing results, not a fresh pipeline reproduction, remote repository audit, or email search. No messages were sent and no analysis/manuscript/package files were changed.

The analysis is substantially prepared, but the September 21 handoff is not complete. Current calibration outputs date to August 7; fossil outputs date to August 6; local manuscript files date to May 26. The working tree was clean at inspection and the latest local analysis commit is `7c9f668` (August 7). Dates are local evidence, not proof that nothing has happened elsewhere.

## Commitments from correspondence

- By September 21: updated calibration validation (Dana calls it Table 1), fossil climate application (Dana calls it Table 6), and explanatory text for an October ten-minute conference talk.
- Talk scope: MAT and MAP. LMA is outside the talk scope, not necessarily outside the manuscript.
- Requested comparison: PIP impute, LM site sp+zero impute, and original Peppe 2011 estimates for calibration; closest non-phylogenetic comparison for fossils.
- Dana supports using each occurrence’s site age, including fossils only identified to genus/family/order. Dana objects to losing within-species climatic variation through cross-site trait averaging and requests explicit disclosure.
- James promised to ensure implementation in dilp and offered to contact Brian. Dana’s thanks acknowledges the offer but does not establish that Brian was contacted. Michelle’s involvement is also unresolved in the supplied thread. No email access was used to verify either.
- Dana did not explicitly settle whether morphotypes separated by long time intervals should be treated as one taxon.

## Calibration results available now

Source: `tables/loso_cv_rmse.csv`.

| Model | MAT RMSE (°C) | ln(MAP) RMSE | Evaluation |
| --- | ---: | ---: | --- |
| PIP impute | 3.414 | 0.525 | 10-fold site-grouped CV; 93 sites |
| LM site sp+zero impute | 4.217 | 0.672 | Same CV; 93 sites |
| Published DiLP coefficients | 3.791 | 0.554 | In-sample reference; 93 MAT / 91 MAP sites |

The published-coefficient row is a current score using fixed published coefficients, not a newly cross-validated competitor or necessarily the exact error statistic printed in Peppe 2011. The comparison must retain that distinction.

LM site sp+zero is a useful practical baseline, but PIP versus this LM changes both fitting/aggregation and phylogenetic treatment. PGLS versus PIP isolates the prediction-time correction; LM species provides another complementary control. The best site LM in these outputs is untoothed-excluded impute (3.722 / 0.583), and the lowest MAP PIP RMSE is complete-case training (0.508 versus impute 0.525). Dana’s simple impute comparison is therefore reasonable but not a comparison exclusively against the best alternatives.

`04_fossil_predictions.R` uses the backward-compatible specimen-impute site model. Inspection of stored `site_models.rds` confirmed its MAT coefficients and imputer match sp_zero_impute; current CV likewise reports identical specimen/sp_zero performance. For a final release, select/name the requested config explicitly instead of relying on this alias.

## Fossil status and implementation gaps

- Main `code/04_fossil_predictions.R:301` grafts each occurrence at its own age, and lines 434–458 use within-site traits plus occurrence-specific phylogenetic adjustments. The age logic applies across genus/family/order/root placements. The pooled fossil species comparator described in the original audit was removed on 2026-09-11; all retained fossil predictions use local traits and occurrence ages.
- Extant PIP training still averages each species across training sites. Held-out and fossil predictions preserve within-site species trait means. Explain this distinction and the loss of within-species information during fitting. Cross-site averaging is a choice of this implementation, not a universal requirement of phylogenetic modelling.
- Two primary occurrence placements failed with `'position' is larger than the branch length`: Cornus nebrascensis at Williston Basin I and Crataegus sp. (rp42) at Republic. The code warns and drops them (`04:322–333`). Thus primary fossil predictions contain 359 of 361 input occurrences. Resolve these before finalizing numbers and counts; the values below are provisional stored results.
- Placement also needs an order-dependence check: `04:165` derives matching genera from every tip in the growing tree, including previously grafted occurrence labels. Later fossils can therefore use earlier fossils as placement evidence. Restricting taxonomic anchors to the intended scaffold and checking shuffled input order should be part of the placement repair; the numerical effect has not been tested in this audit.
- Primary formal-only placement logs contain 191 successful root placements. These are distinct from failed placements. Taxonomic resolution and formal-only versus provisional quoted taxonomy should be documented.
- `release/code/04_fossil_predictions.R` still uses the older species-mean-age placement for the site model.
- Local `dilp/R/dilp_pgls.R:162–203` also averages ages before placement and reuses species corrections at lines 217–244. It preserves within-site prediction traits but does not use occurrence-specific tree tips. It lacks the main pipeline’s fossil bagged imputation and formal-only quoted-taxonomy handling. The local package is clean at pinned commit `b29be90` (April 30); the staged update files date to May 4. Updating constants alone will not synchronize behaviour.
- Actual input has 347 species labels across 361 occurrences. Twelve labels recur across sites, eleven across different ages, with maximum age span 8.5 Ma. The email’s Fox Hills–Republic ~17 Ma scenario is not present as a repeated species label in this input. Taxonomic validity remains a collaborator question.

## Provisional fossil comparison

Source: `tables/fossil_site_comparison.csv`, formal-only, impute variants. MAP is cm. These are existing outputs, not new calculations of predictions.

| Site | MAT PIP | MAT LM site | MAP PIP | MAP LM site |
| --- | ---: | ---: | ---: | ---: |
| Fox Hills | 16.3 | 18.1 | 156 | 83 |
| Williston Basin I | 17.7 | 14.6 | 174 | 138 |
| Williston Basin II | 16.8 | 13.8 | 166 | 102 |
| Palacio de los Loros PL1 | 17.1 | 12.9 | 170 | 118 |
| Palacio de los Loros PL2 | 17.6 | 17.6 | 169 | 83 |
| Williston Basin III | 16.6 | 16.1 | 162 | 106 |
| Cerrejon | 21.3 | 25.9 | 223 | 288 |
| Hubble Bubble | 19.2 | 17.8 | 170 | 96 |
| Laguna del Hunco | 16.7 | 11 | 169 | 139 |
| Republic | 14.3 | 8.6 | 150 | 67 |
| Bonanza | 17 | 10.3 | 151 | 116 |

## Remaining manuscript work

- Refresh Table 1 and derived prose in `doc/model_validation.md:55`; fix 92 versus 93 sites, outdated complete-case sample descriptions, and outdated imputation language. The current code refits the species imputer inside each CV fold (`03_loso_cv.R:365`), contrary to the document’s claim that it is fitted globally before CV.
- Rename the old “Peppe” fitted LM configuration to “untoothed excluded,” distinguishing it from actual published DiLP coefficients.
- Replace stale fossil values/counts in `doc/empirical.md:23` after resolving placement failures, and add the LM comparison. The local fossil climate table is numbered Table 4, whereas Dana refers to Table 6: reconcile against the assembled manuscript before renumbering.
- Remove the outdated claim that every occurrence of a species receives the same phylogenetic adjustment (`doc/empirical.md:15`, `doc/fossil_predict.md:16`, main README). Explain training grand means versus occurrence-specific prediction traits and ages.
- Document formal-only primary taxonomy and informal sensitivity. Stored rounded sensitivity differences reach 0.8 °C and 5 cm MAP for PIP site.
- Additional consistency check before interpretation: CV averages species predictions on the log(MAP) scale (`03:746`), whereas fossil predictions exponentiate before averaging (`04:455–458`). State or reconcile the site-level estimand before using the CV error as direct uncertainty for the fossil arithmetic-mean estimate.
- Regenerate Word versions after Markdown changes and prepare a concise MAT/MAP-only talk extract. Keep LMA work outside this handoff unless separately requested.

## Suggested order to meet September 21

1. Resolve the two occurrence graft failures and verify primary site ages, retention counts, and the explicit LM baseline.
2. Synchronize occurrence placement, taxonomy and missing-trait behaviour across main analysis, release, and dilp; refresh package model objects as needed and verify prediction parity.
3. Regenerate affected fossil outputs; rerun calibration only if training/CV changes warrant it. Reconcile the MAP aggregation interpretation.
4. Update manuscript tables/methods and produce Dana’s short comparison, clearly separating published in-sample scores from CV.
5. Send the completed update to Dana by September 21. Brian outreach is settled per James's September 7 follow-up. Sending remains a separate action requiring James’s instruction.

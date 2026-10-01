# Decisions

## 2026-09-25: stage rewrite scope

| Question | Decision | Reason |
|---|---|---|
| What "unbiased" means | Same tests, thresholds, plots and tables at every level and in both arms. Fixes to methods that favour one arm are deferred. | Nothing is chosen by hand; every nominal hit sits beside the count chance predicts. |
| What "comprehensive" means | Protein, pathway and module each run the nine contrasts, the eight classification tasks and phenotype association, plus their own extras. | Parity across levels, keeping filtering QC, fgsea collapse, volcano rings, WGCNA preservation and STRING. |
| Phenotype link | A sample-level limma model per trait (abundance ~ timepoint + trait, subject blocked, trait split into within- and between-person terms), plus the change-score correlations. | One model answers both questions; the change score is a second view of the within-person question. |
| HR/LR pair imbalance | Deferred. HR has 6 training and 7 acute pairs, LR 8 of each. | Every caption and sub-stage README states n per task, the smallest attainable p and the chance count. |
| Step layout | Every level runs contrasts, then classify, then associate, then its extras. | One level teaches the others. |
| PDF scope | Every nominal hit gets a panel, 12 to a page, by p; chance line on the first page. | Complete, with no hand-set cutoff. |
| Stages in scope | 00 to 04. `05_Summary` is a short plan in its README. | 00 builds the per-biopsy phenotype table; 01 adopts the output conventions. |
| Repeated logic | Each script carries its own short copy, written the same way at every level. | Scripts stay self-contained and lean on package calls. |
| Script format | Plain R scripts, run with `Rscript`. | One run command, no Quarto, no ignored HTML. |
| Checks | Snapshot every workbook, rerun, diff; report each number that moves and why. | Nothing moves without a stated reason. |
| Per-biopsy traits | fCSA mixed, type I and type II; mCSA; 1RM leg press and extension; fibre counts mixed and type I; % type I fibres. Percent change per trait joins the change-score outcomes. | % type I carries the positive control; percent change is one value per person. |

## 2026-09-25: added during the build

| Question | Decision | Reason |
|---|---|---|
| Plot code | Package wrappers only: `pROC::ggroc()` for ROC pages, `ggpubr::ggscatter()` for change-score scatters, plain `ggplot2` elsewhere. Pages of twelve come from `split()`, not `ggforce::facet_wrap_paginate()`. | Paginate rebuilds every panel on every page; on 400-page PDFs that is too slow. |
| Inputs | `00_Input` holds hand-curated study data only: `metadata.csv` (one row per MS sample) and `phenotype.csv` (one row per subject, pre and post side by side). Builders and derived inputs are gone. | Same pattern as BFR_limpa. Change scores from the new file equal the old ones to 1e-9. |
| Observation floor in the sample model | A protein enters a trait's model only when two thirds of that trait's biopsies saw it. | TIGAR, seen in 15 of 30 biopsies, reached BH < 0.05 on between-person mCSA through half-missing data; the change score already used this floor. |
| Derived trait name | `type1_share`, the type I fibre count over the mixed count, in percent. | `pct_type1` would have produced `pct_pct_type1` for its percent change. |

## 2026-09-27: nominal reporting audit

| Question | Decision | Reason |
|---|---|---|
| Chance ratio | Removed everywhere. Tables report tested, p < 0.05 and BH < 0.05 counts; the smallest attainable p stays. | 5% of tests is not the null count here: the Wilcoxon p is discrete with 6 to 8 pairs, and proteins and sets are correlated. This supersedes the chance entries above. |
| Paired tasks | Within-arm tasks are drawn as paired grids with the signed-rank p and median paired change; ROC grids only for the four HR-vs-LR tasks. | AUC and ROC treat paired biopsies as independent groups, so they described a different comparison from the signed-rank p. |
| Set contrast pages | Every fgsea collapse survivor at FDR < 0.05 and every fry set at p < 0.05, not the ten strongest. | A top-ten cut is a hand-set cutoff. |
| Summary figures | The scratch summary figures are dropped; `05_Summary` stays a plan. | They carried chance dots and hand-picked examples. |

## 2026-09-27: plot pass

| Question | Decision | Reason |
|---|---|---|
| Baseline hits | Drawn like every other task; the task label drops "(floor)". The contrast keeps its floor role in `design.rds`. | Every nominal hit gets a panel; leaving one task out was a hidden cut. |
| Change-score scatter | ggplot2 `geom_point()` with a least-squares line per arm, captioned as a guide; the tests stay Spearman. `ggpubr` leaves the scripts. | `ggscatter()` did nothing ggplot2 does not, and 25 packages entered the lock only through it. This supersedes the plot-code entry of 2026-09-25 for scatters. |
| Overviews | Each contrasts, classify and associate PDF opens with ComplexHeatmap pages of every feature: biopsies at the contrasts step, effects beside a status map (not significant, p < 0.05, BH < 0.05, not tested) at every step. Pathway pages repeat per collection. | The grids show hits one at a time; the overviews show every feature at once. ComplexHeatmap was already in the lock through proteoDA. |
| Overview rows | Proteins by average-linkage clustering on 1 - Pearson r of the imputed matrix, sets the same on singscore within collection, modules by size. Display only. | One order per level, the same in every step. |
| Panel order | Within each comparison: p, then the larger effect, then feature ID. | For a fixed n, p is monotone in the effect, so the strongest come first and ties are deterministic. |
| Panel n | Each header gives the feature's own n: complete pairs, HR and LR subjects, biopsies and subjects, or subjects. | Missing values change n per feature. |
| Plot index | Every classify and associate workbook gains `plot_index`: result key, PDF page, panel position. | One lookup from a table row to its panel. |

## 2026-09-29: responder analysis review

| Question | Decision | Reason |
|---|---|---|
| Training interaction | Keep it primary; add an ANCOVA (`T2 ~ T1 + arm + blood`, and `T3 ~ T2 + arm + blood`) at the protein, set (singscore) and module levels. | Every interaction hit shows HR starting low and catching up, or starting high and dropping (OXPHOS NES −2.59 at baseline, +1.82 interaction); ANCOVA separates that from regression to the mean. |
| Effect size | Nominal p < 0.05 with logFC reported; no `treat()` or fold filter for now. | At this n a nominal hit is already a large observed effect (smallest Training_Interaction hit 1.35-fold). |
| Counting agreement | Rank candidates by raw nominal-hit count; `read_me` notes the training-axis tests share data. | User's choice. |
| Traits | comp_hypertrophy, fCSA mixed, type I, type II, mCSA, leg-press and leg-extension 1RM. Fibre counts and `type1_share` leave the associate steps. | Size and strength outcomes only. A fibre-type pilot found no HR/LR difference; it was deleted. |
| Change unit | Percent change only; the `d_*` absolute screens go. | Comparable across traits; one screen per trait. |
| Per-trait splits | None. HR/LR is the only grouping; every trait runs continuously. | fCSA splits match HR/LR in 14 of 16 subjects; leg press in 8 of 16. |
| Strength | Two-predictor model per level: training change ~ % comp_hypertrophy + % strength gain, once for leg press, once for leg extension. | Strength and hypertrophy are unrelated here (4/4/4/4 crosstab); leg-press and leg-extension gains correlate −0.45, so no composite. |
| Blood | Loading-adjusted blood index (Hb minus sample median log2, centred) as one term in every model; the unadjusted fit is kept as a sheet. | LR T3 carries +1.25 vs about 0 elsewhere. Pilot: training contrasts move r ≥ 0.94; Acute_LR OXPHOS goes from −1.08 to +2.00*; the acute interaction sharpens. |
| Set tests | Report fgsea and fry side by side; the README explains fry (samples as the unit, direction, power set by people). | User's choice. |
| Quantification | Stay on proteoDA; mirror BFR's layout, not limpa. | Only a protein-level matrix exists; no DIA-NN precursor report. |
| Naming | Keep V2 step names; add a "maps to BFR" table to the README. | One slot layout at every level. |
| Summary | `05_Summary` builds a candidate table only, no figures. | User's choice. |

Open: whether `volume_load_total_kg` stays as its own screen; whether the blood, ANCOVA and two-predictor helpers move into a tested `R/` file.

## 2026-09-29: protein candidate rules

| Question | Decision | Reason |
|---|---|---|
| Organisation | Two lists, training and acute; a protein on both is flagged. | Chronic remodelling and the bout response are different questions with different caveats. |
| Training-list DE evidence | Training_Interaction, ANCOVA training, Trained_HRvLR, Acute_Interaction. | User's choice; Acute_Interaction counts on both lists. |
| Acute-list DE evidence | Acute_Interaction, ANCOVA acute, Acute_HRvLR. | |
| DE hit | Nominal p < 0.05, logFC reported. | FDR is empty on every HR-vs-LR contrast even at 0.25. |
| Classification hit | Wilcoxon p < 0.05, AUC and median change reported. Training list: training change and T2; acute list: acute change. | At 6 to 8 per arm p < 0.05 already needs AUC of about 0.83. |
| Association | Pooled, within HR, within LR, and the trait-by-arm slope difference, over the seven traits' percent change; plus the two-predictor growth and strength terms. | Within-arm links and a formal test that they differ. |
| Counting | Every p < 0.05 counts once. | Simple and fully reported; correlated traits are not grouped. |
| Inclusion | Proteins with at least one hit, sorted by hit count then smallest p. Baseline logFC and p, blood-unadjusted p and module are columns, never filters. | Nothing hidden. |

## 2026-09-29: how the review was built

| Question | Decision | Reason |
|---|---|---|
| Blood covariate | Each protein's slope on the loading-adjusted blood index, fitted beside the six group means with subject blocked, comes off in `02_normalize.R`; `blood` stays in the design for its degree of freedom. | One adjustment reaches every method (limma, fry, singscore, WGCNA, rank tests) without a covariate in each. Refit counts match the covariate model contrast for contrast; the `unadjusted` sheet reproduces the old fit exactly. |
| Modules | WGCNA on the adjusted matrix finds 15 modules, not 12. | Blood co-variation no longer binds proteins together. |
| Leg-press slope difference | Kept as reported: 27 proteins at FDR < 0.05, 26 in the contractile module (pink). All 27 stay nominal without HR_S03 (+108% leg press) and with the outcome ranked; 15 stay at FDR < 0.05 without HR_S03. | Checked for leverage before reporting. |
| Classification wording | Per-feature Wilcoxon and AUC is single-feature discrimination, not a predictive classifier; no AUC is out-of-sample. | Literature check: Bamber 1975, Hanley & McNeil 1982 (AUC equals Mann-Whitney), Christin et al. 2013, Ambroise & McLachlan 2002, Varma & Simon 2006. |

## 2026-09-30: panel, extra tests and shared code

| Question | Decision | Reason |
|---|---|---|
| AUC intervals | DeLong 95% interval on between-arm tasks, blank at AUC 0 or 1. | DeLong's variance is zero at perfect separation and the interval collapses to a point; at this n it runs narrow. |
| Panel inputs | Five panels: T1, T2, T3, training change, acute change. | T1 is the only one that predicts in advance; the rest describe. |
| Panel features | Module eigengenes and all proteins, run separately. | A panel that works only on 1,900 proteins is a warning of overfitting. |
| Panel method | `nestedcv`: elastic net, alpha 0.5 fixed, `lambda.1se` inside each fold, leave-one-out outer loop, imputation and scaling inside folds, 1,000 whole-pipeline label permutations, p = (b + 1)/(m + 1), selection frequency per feature, labelled exploratory with no external validation. | Airola et al. 2011; Smith et al. 2014; Varma & Simon 2006; Cawley & Talbot 2010; Meinshausen & Bühlmann 2010; Phipson & Smyth 2010; TRIPOD+AI 2024. Leave-pair-out would be less biased but no standard wrapper offers it. |
| Eigengene panel | Modules come from all 16 subjects; the leak is stated in the output, not removed. | Rebuilding WGCNA inside every fold is impractical. |
| Extra checks | Leave-one-subject-out on the key tests; trace why Training_LR fgsea sets fell from 121 to 4 after blood adjustment. | Single-subject results; largest shift on the training side. |
| Code structure | One `R/helpers.R` sourced by every script, unit-tested with testthat; row-wise packages (`matrixTests`, `limma`) replace hand loops. Supersedes "each script carries its own copy" (2026-09-25). | One fix lands everywhere; tests become possible. |
| Build order | Refactor first with no number moving (workbook diff), then the new checks, then the panel. | New work builds on the shared helpers. |

## 2026-09-30: Training_LR after blood adjustment

| Question | Finding | Reason |
|---|---|---|
| Why Training_LR fgsea sets fell from 121 to 4 | Superseded 2026-10-01: the fall came largely from the HBG2-inflated blood index. Original reading: expected. LR's blood index falls 0.63 units from T1 to T2, HR's 0.21, so the unadjusted contrast carried less blood after training. Protein logFC barely moves (r = 0.94, median shift 0.05); set NES agree at r = 0.84. The lost sets are 88 Reactome sets, mostly one proteasome-driven cluster (APC/C, neddylation, CLEC7A, BCR signalling, DNA damage checkpoints), whose NES shrank from about -2.1 to -1.8 and crossed FDR 0.05 together. Proteasome subunits carry a larger blood slope than the proteome (median 0.021 against 0.002); red cells hold proteasome. | Threshold effect on a correlated cluster, not a failure of the adjustment. |
| Refactor check | 21 workbooks rerun after moving shared code to `R/helpers.R`: every result sheet identical; only `.xlsx` input checksums differ, because writexl stamps a creation time. | No number moved. |

## 2026-09-30: panel outer loop

| Question | Decision | Reason |
|---|---|---|
| Outer cross-validation | Stratified 5-fold, 10 repeats, replacing leave-one-out. | On T1 proteins `lambda.1se` kept nothing in every fold, and pooled leave-one-out then scored AUC 0.00: each held-out HR subject leaves a training set tipped towards LR (Airola et al. 2011). Replicated CV is unbiased in Smith et al. 2014. |
| Reading the AUC | Against the permutation null only. | A 5-permutation smoke run put the null median near 0.39; with 7 against 8 no fold split is balanced. |
| Missing values | Proteins observed in every subject of the input; nothing imputed. | glmnet takes no NA, and imputation outside the folds would leak. 974 to 1,115 proteins per input. |

## 2026-09-30: panel results and the T3 blood confound

| Question | Finding | Reason |
|---|---|---|
| Panels beyond their nulls | Modules at T3: AUC 0.83 against a null median of 0.38, permutation p = 0.006, carried by salmon (kept in 100% of models) and yellow (96%). Proteins in the acute change: AUC 0.58, p = 0.042, RPS20 kept in 64%. The other eight panels sit inside their nulls. Ten panels were tested, so p = 0.042 is within chance and p = 0.006 is near a Bonferroni line of 0.005. | 1,000 permutations each. |
| Reading the T3 module panel | Confounded with blood. At T3 the salmon eigengene correlates rho = -0.83 with the blood index (p = 1e-4) and yellow -0.63; both sit lower in LR, the bloodier arm. The blood adjustment fits its slope within group cells, so the HR-against-LR blood difference at T3 is aliased with the group means and cannot be removed. | Any HR-against-LR result at T3 carries this caveat, not only the panel. |

## 2026-09-30: exploratory corrections, tests removed, outputs audited

| Question | Decision | Reason |
|---|---|---|
| Correction beside nominal results | BH and Storey's q in every table, within the same screen; the candidate tables carry both per hit and count hits below 0.20. | User's choice. Storey uses a bootstrap pi0, which copes with discrete Wilcoxon p; screens under 100 tests (modules) set pi0 = 1, so q equals BH. |
| Unit tests | `tests/` removed at the user's request. Supersedes the tested-helpers part of 2026-09-30. | Overrides the global testing policy for this project, by explicit instruction. |
| Output audit | All 58 tracked files under `c_data/` and `b_reports/` are written by their own step's script; two empty `b_reports/` folders left by older versions were removed. | |

## 2026-10-01: statistical audit and limpa rebuild

Three read-only audits (preprocessing; protein and pathway statistics; network and summary) checked
every stage against its package documentation and its outputs. The interview settled each finding.

| Question | Decision | Reason |
|---|---|---|
| Blood index | Mean log2 of HBB and HBA1 against the sample median; HBG2, HBD and HBG1 dropped. | HBG2 was missing in 20 of 48 runs, the low-blood ones, and averaging what was present lifted those runs by about 1.57 log2. With the fixed index the arm-by-T3 blood difference falls from 1.21 (p = 0.023) to 0.84 (p = 0.23). |
| Blood handling | A covariate in every limma model, no pre-subtraction; `removeBatchEffect` for methods that cannot take one (group design; time point only for the panel, so no label leak). | Standard covariate route (limma user guide); the earlier pre-subtraction used unshrunk slopes. BFR_limpa uses sample weights instead. |
| Missing values | limpa row-wise: `dpc()`, `dpcQuantByRow()`, `dpcDE(block = subject)`. Replaces the proteoDA fit and missForest. | Missingness is intensity-dependent (rho -0.77); 17 of 49 primary hits rested on a cell with 3 or fewer values; missForest imputed mid-range and shrank variance. The limpa vignette gives this route for protein-level input; Smyth (support.bioconductor.org/p/9163224) requires intensity-based log values. DPC slope 0.67. |
| Set tests | fgsea primary; cameraPR added as a corroboration column with a `camera_agrees` flag; fry kept as the strict test. | cameraPR kept 92% or more of fgsea hits nominal in the same direction (rank rho 0.94-0.97), so it trims without washing out; it also finds sets on Baseline_HRvLR, so that contrast is a comparison, not a null. |
| ANCOVA | Blood-corrected observed values, variances moderated with `squeezeVar`. | Residual df were 3 to 12. |
| Classification | Unchanged. | User's choice. |
| Arm slopes | Five observed subjects per arm; `comp_hypertrophy` out of the arm-difference screen. | The overall floor let an arm's slope rest on 3 or 4 people; arm is the median split of comp_hypertrophy. |
| Modules | Zsummary bands per Langfelder 2011 (2-10 weak to moderate); random-halves calibration; ICC with time point; kME on the subject-centred matrix. | The calibration splits all 16 subjects into random halves of 8 rather than each arm in half: same network size as the arm comparison, one population. |
| Candidates | Raw hit-count ranking kept; ANCOVA label bug fixed; cameraPR not counted. | User's choice; the bug left every ANCOVA hit without a stability call. |
| Chance line | `05_Summary/01_Permutation`: 200 shuffles of subjects' phenotype rows at all three levels, context only. | The baseline list is not a null. |
| Step numbering | `05_Summary` reordered to run order: 01_Permutation, 02_Candidates, 03_Panel. | Candidates read the permutation scores. |

## 2026-10-01: results after the limpa rebuild

| Measure | Before | After | Reading |
|---|---|---|---|
| Training_Interaction nominal proteins (median abs t) | 49 (0.63) | 41 (0.58) | Every protein now tested in every contrast; calibration unchanged or tighter. |
| Acute_Interaction nominal (median abs t) | 99 (0.72) | 79 (0.69) | |
| Interaction hits stable to leaving one subject out | 18 training, 23 acute | 5 training, 11 acute | limpa's weights lean less on a few sparse values. |
| ANCOVA nominal (training, acute) | 62, 173 | 53, 171 | Now moderated. |
| fgsea FDR < 0.05, Training_LR | 4 | 116 | The earlier collapse (entry of 2026-09-30) came largely from the HBG2-inflated blood index over-correcting, not from proteasome biology. |
| fgsea hits cameraPR corroborates, Training_Interaction / Acute_Interaction | n/a | 62 of 64 / 189 of 206 | |
| Modules | 15 | 10 | Built on limpa's blood-corrected matrix. |
| Arm networks inside the random-halves preservation range | n/a | 35 of 37 | HR and LR share structure as fully as one population does with itself. |
| Leg-press slope difference, FDR < 0.05 | 27 | 20 | Per-arm floor of five. |
| Candidate lists against the permutation chance line | n/a | every list at chance (share of shuffles at or above the real count 0.10 to 1.0) | A random relabelling yields as many multi-hit candidates as the real labels. |
| Panels | modules T3 p = 0.006 | all at chance (p 0.28 to 0.99), no feature kept in half the models | The T3 result was carried by blood and the label-informed correction. |

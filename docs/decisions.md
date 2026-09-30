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

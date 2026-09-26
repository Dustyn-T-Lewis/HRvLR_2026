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

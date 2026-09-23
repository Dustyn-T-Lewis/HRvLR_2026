# HRvLR_V2

Skeletal-muscle proteome of 16 subjects, 8 labelled High and 8 Low Responder by
a median split of a composite hypertrophy score, biopsied at baseline (T1),
after training (T2) and after an acute bout (T3). This project asks whether the
proteome separates the two arms, and whether it tracks how much each subject
adapted, at three feature levels: proteins, pathways and co-expression modules.

## Layout

```
00_input/          raw matrix, sample metadata, phenotype table, GO-Slim .obo
01_Preprocess/     proteoDA filtering, cycloess normalization, missForest
02_Proteins/       the nine contrasts, protein screens, protein packet
03_Pathways/       gene sets and themes, fry and fgsea, singscore, screens, packet
04_Networks/       WGCNA modules, their characterisation, screens, packet
05_Figures/        planned; reads only c_data from 02-04
functions/         shared code, flat
run_all.R          rebuilds 01 to 04 in order
```

Each stage builds its matrices first, runs its tests on them, and ends in a
captioned PDF packet. Every unit carries `a_script/` `b_reports/` `c_data/`.
`c_data/` is tracked and holds one handoff `.rds` and one `.xlsx`; renders in
`b_reports/` are not. Paths go through `here::here()` and every stochastic step
is seeded with 42.

## Pipeline

**01_Preprocess.** proteoDA end to end, as in YvO: `DAList` →
`zero_to_missing` → a four-method outlier consensus (drops S29_T1, S28_T2,
S28_T3) → `filter_proteins_by_group` → `filter_samples` → `write_norm_report`
→ `normalize_data(norm_method = "cycloess")`. Contaminants are removed by
protein identity (HPA tissue, a curated blood list). The result is 1900
proteins × 45 samples. missForest fills that matrix for the steps that need it
complete (fry, singscore, WGCNA); nothing else reads the imputed values.

**02_Proteins.** `add_design("~ 0 + group + (1 | subject)")`,
`add_contrasts()`, `fit_limma_model()`, `extract_DA_results()`,
`write_limma_plots()`. Nine contrasts: four within-arm changes (training
T2−T1 and acute T3−T2, per arm), the arm difference at each timepoint, and two
interactions. The fit reproduces V1's committed numbers to 5.6e-12
(`02_verify_v1.R`). The handoff `proteins.rds` carries protein × contrast
matrices of logFC, t, p and BH, the abundance matrix, and three subject
windows (T1 level, training change, acute change).

**03_Pathways.** MSigDB 2026.1 Hallmark, Reactome and GO:BP plus one set per
GO-Slim generic term, rewritten to protein ids and kept at 15 to 500 detected
members: 1596 sets. Each GO:BP set is themed by its most specific GO-Slim
ancestor (657 of 1156 fall under one). `limma::fry` tests every set in every
contrast on the protein design and subject block; fgsea on the moderated t
supplies NES for display, with `collapsePathways` marking non-redundant sets.
singscore gives the set × sample matrix the screens read.

**04_Networks.** WGCNA, signed, bicor, defined on abundance centred within
subject and scored on raw abundance (`functions/shared_wgcna.R`): power 8, 12
modules, 265 proteins unassigned. Each module is characterised by
clusterProfiler ORA against the 03 sets (universe the 1900 detected proteins),
its ten highest-kME hubs, its dominant GO-Slim theme, and its STRING v12
edge density against a null that shuffles module labels. Eigengenes are then
fitted on the nine contrasts with the protein design.

**The screens** (`functions/classify.R`, the same code at all three levels):

- *Classification.* Per-feature AUC (pROC) and Wilcoxon p on seven tasks that
  mirror the contrasts: T2 vs T1 and T3 vs T2 within each arm (paired), and HR
  vs LR on the T1 level, the training change and the acute change.
- *Association.* limma with a continuous predictor, one row per subject, on
  the three windows × ten phenotypes.

BH runs within each task, window × phenotype cell, or collection, never across
them. Every table reports the nominal count beside the count chance predicts.
No sweep-level correction is applied; results are exploratory.

## Current results

| Level | Contrasts, BH < 0.05 | Classification, BH < 0.05 | Association, BH < 0.05 |
|---|---|---|---|
| Proteins (1900) | 0 of 9 contrasts; lowest q 0.071 (Acute_HR) | 0 | 0 |
| Pathways (1596, fry) | 8 set-contrast pairs, all acute | 0 | 1 |
| Modules (12) | 0 | 0 | 1 |

- **Proteins.** 959 nominal and 464 pi-score calls across the nine contrasts,
  against roughly 95 nominal per contrast expected by chance. The
  within-subject correlation is 0.189.
- **Pathways.** fry calls five sets in Acute_HR (Hallmark G2M checkpoint and
  MYC targets up, fatty-acid metabolism and adipogenesis down, GO-Slim mRNA
  metabolism up) and Hallmark heme metabolism in Acute_LR, Acute_HRvLR and
  Acute_Interaction. The one association is Reactome basal-body anchoring
  against the acute change and `comp_hypertrophy`.
- **Modules.** Six modules share more STRING edges than the null (BH < 0.05):
  magenta (translation, 16×), tan (respiratory chain, 19×), pink (striated
  muscle contraction, 11×), purple (translation initiation), yellow (aerobic
  respiration) and brown (cytoskeleton). The one association is greenyellow,
  the extracellular-matrix module, against the acute change and `d_mcsa`
  (BH 0.037 across 12 modules). An ECM module against whole-muscle CSA was
  also the closest result in the earlier continuous design; see History.

Twelve modules make BH a much weaker filter than 1900 proteins do. Under a
global null each of the 30 association cells has about a 5% chance of at least
one BH hit, so about 1.5 cells with a hit are expected by chance; one was
observed at the module level and one at the pathway level. Each packet states
its own chance line.

## Phenotypes

Ten, all T2−T1 changes except `volume_load`: `comp_hypertrophy`, fibre CSA
(`d_fcsa_I`, `d_fcsa_II`, `d_fcsa_mixed`), MyoVision fibre counts
(`d_nfibre_mixed`, `d_nfibre_I`, counts despite "fCSA" in their source names),
whole-muscle CSA (`d_mcsa`), 1RM (`d_1rm_legpress`, `d_1rm_ext`), and total
kilograms lifted (`volume_load`). They form roughly five independent axes. The
HR/LR label is the exact median split of `comp_hypertrophy`, so it separates
that composite's fibre-CSA ingredients by construction.

## Running it

```sh
Rscript setup.R      # restore the renv library once
Rscript run_all.R    # 01 to 04, one R session per step, about 3 minutes
```

`run_all.R` logs each step to `.runlogs/`. The packets land in
`02_Proteins/03_Packet/b_reports/`, `03_Pathways/05_Packet/b_reports/` and
`04_Networks/04_Packet/b_reports/`. Tests: `Rscript tests/testthat.R`.

STRING v12 (`9606.protein.links` and `.aliases`) is read from
`00_input/downloads/`, which git ignores; download both from string-db.org
before running 04.

## History

V1 tested HR versus LR at four levels and found nothing. V2 has since tried
three framings before this one, each recorded in git history with its code and
figures in the untracked `archive/`:

- **Classification** (`f4de361` to `0fd4175`): ten candidate labels and a blind
  subtype search. The label separates only its own ingredients, and no
  clustering beat a no-cluster null.
- **Continuous** (`5e9211b` to `158ec3d`): 3 levels × 6 windows × 10
  phenotypes. No cell survived correction across the sweep; the closest result
  was an ECM module against change in whole-muscle CSA.
- **Tertiles** (`f61e280`): three response groups. This pass also caught
  sparse proteins producing false hits, which the missingness filter now
  handles.

The current spine restarts on the nine contrasts (`6c93f0d`) and adds the
classification and association screens at every level.

## References

Atkinson G, Batterham AM (2015). True and false interindividual differences in
the physiological response to an intervention. *Exp Physiol* 100:577-588.

Xiao Y, Hsiao TH, Suresh U, et al. (2014). A novel significance score for gene
selection and ranking. *Bioinformatics* 30:801-807.

# 05_test — GSVA pathway scores + mixed models (experimental)

Exploratory stage testing whether single-sample pathway activity, modelled with
subject random effects, adds to the contrast-based concordance work in `04_`.
**Not a manuscript figure yet** — a sandbox to judge the method.

## Pipeline

1. `setup.R` — protein matrix → gene-symbol expression (`limma::avereps`), Hallmark
   + GO Slim gene sets, design metadata, response phenotypes.
2. `01_gsva_scores.R` — per-sample pathway activity via GSVA (Gaussian kcdf for
   log-intensity data). Scores computed once on the full matrix.
3. `02_lmm_trajectory.R` — per pathway `score ~ group * timepoint + (1|subject)`;
   the group × timepoint term is the HR-vs-LR trajectory divergence. BH across
   pathways. Figures: trajectories of the top pathways, interaction ranking.
4. `03_lmm_phenotype.R` — `score ~ phenotype * timepoint + (1|subject)` (does the
   trajectory track the response), and baseline (T1) score → phenotype with
   leave-one-out CV referenced to a permutation null.
5. `04_multidb.R` — repeat GSVA and ssGSEA across all five databases, fit the
   group × timepoint LMM per pathway per method, and correct within each database
   (each has its own size and redundancy), keeping nominal p alongside the FDR.
6. `05_parallel_figures.R` — method-agreement and per-database summary figures
   comparing GSVA and ssGSEA.

`run.R` sources `setup.R` then `01`–`05`.

## Method validity (verified citations)

- GSVA: Hänzelmann, Castelo, Guinney 2013, BMC Bioinformatics 14:7 (PMID 23323831).
- ssGSEA: Barbie et al. 2009, Nature 462:108 (DOI 10.1038/nature08460);
  GSEA: Subramanian et al. 2005, PNAS 102:15545 (PMID 16199517).
- lme4: Bates et al. 2015, J Stat Softw 67(1); lmerTest: Kuznetsova et al. 2017,
  J Stat Softw 82(13) — Satterthwaite df at small n.
- LMM for longitudinal / per-feature omics: Mallick et al. 2021 (MaAsLin 2),
  PLoS Comput Biol 17:e1009442 (PMID 34784344); Taheriyoun et al. 2026,
  Comput Struct Biotechnol J 31:301.
- FDR: Benjamini & Hochberg 1995, JRSS-B 57:289.
- CV / selection bias: Ambroise & McLachlan 2002, PNAS 99:6562.

## Caveats (do not skip)

- **GSVA scores are cohort-relative** (Hänzelmann 2013; Bhuva et al. 2020, NAR
  48:e113) — they do not transport to new samples unchanged, which limits any
  out-of-cohort prediction claim.
- **Prediction is the weak point**: n=16; a strictly honest LOO would re-score
  GSVA inside each fold, and the best-of-many-pathways Q² is selection-biased.
  The permutation null is the reference that keeps the claim honest.
- No published precedent for GSVA + LMM in muscle/exercise omics — defensible as a
  composition of validated parts, framed as exploratory.

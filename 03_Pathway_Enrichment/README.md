# 03 · Pathway Enrichment

Tests each gene set on each contrast, scores every sample on every set, and tests the scores.

| Step | Runs | Writes |
|---|---|---|
| [`00_build_gene_sets`](00_build_gene_sets/README.md) | frozen MSigDB snapshot, protein-to-gene map, size filter, GO Slim sets | `gene_sets.rds` |
| [`01_run_fgsea_and_fry`](01_run_fgsea_and_fry/README.md) | fgsea and fry per contrast, `collapsePathways` | `set_tests.rds` |
| [`02_enrich_volcano_fgsea`](02_enrich_volcano_fgsea/README.md) | protein volcanoes with pathway rings | 8 volcanoes |
| [`03_enrich_scatter_fgsea`](03_enrich_scatter_fgsea/README.md) | HR NES against LR NES, training and acute | 4 composites |
| [`04_run_singscore`](04_run_singscore/README.md) | per-sample set scores | `singscore.rds` |
| [`05_classify_and_associate_sets`](05_classify_and_associate_sets/README.md) | set classification and phenotype association | every nominal set drawn |

## 1,378 sets pass the size filter

The first run fetches four MSigDB collections (release 2026.1.Hs) and writes an RDS with an md5 to
`00_build_gene_sets/c_data/cache/`. Later runs verify the checksum and need no network. Restore the
RDS and its `.md5` together; move both aside to rebuild. `goslim_generic.obo` sits beside it with
its own `.md5`; the script does not fetch it. `00_Input/PROVENANCE.md` has its source and refresh
steps.

A set is tested with 15 to 500 source genes and at least 15 measured. Each GO Slim set holds every
measured gene annotated to the slim term or any GO:BP term below it, and its size rule reads
measured size. All 1,900 proteins carry one distinct symbol, so every protein represents its gene.

| Collection | In source | Tested | Median measured |
|---|---:|---:|---:|
| Hallmark | 50 | 33 | 34 |
| KEGG_Legacy | 186 | 62 | 23.5 |
| Reactome | 1,839 | 328 | 33 |
| GOBP | 7,538 | 897 | 25 |
| GO_Slim | 71 | 58 | 97 |
| Total | 9,684 | 1,378 | |

## fry finds nothing; fgsea calls more sets on the floor than on the primary contrast

fgsea asks whether a set sits at one end of the protein ranking (moderated t, seeded). fry asks
whether the set moved at all under the fitted design, subject block and within-subject correlation
(0.176, estimated on the imputed matrix, since fry takes no missing value). fry tests all 1,378
sets; fgsea tests 1,369 to 1,374 per contrast, because a protein untested in a contrast leaves that
ranking. Each method's own BH runs within each contrast over the five collections pooled.
`collapsePathways` then marks non-redundant fgsea hits in the `main` column and deletes nothing.

| Contrast | fgsea | after collapse | fry |
|---|---:|---:|---:|
| Training_Interaction *(primary)* | 27 | 13 | 0 |
| Acute_Interaction *(secondary)* | 110 | 44 | 0 |
| Training_HR | 68 | 23 | 0 |
| Training_LR | 121 | 18 | 0 |
| Acute_HR | 279 | 83 | 0 |
| Acute_LR | 42 | 18 | 0 |
| Trained_HRvLR | 0 | 0 | 0 |
| Acute_HRvLR | 175 | 58 | 0 |
| Baseline_HRvLR *(floor)* | 60 | 26 | 0 |

fgsea calls 60 sets on the floor and 27 on the primary contrast, so no fgsea list is a finding.

HR and LR NES correlate at rho 0.34 over 1,369 sets for training and 0.35 over 1,370 for the acute
bout.

| Pair | Population | Sets | rho | Significant | Discordant |
|---|---|---:|---:|---:|---:|
| training | all collections | 1,369 | 0.34 | 176 | 15 |
| training | collapse survivors | 38 | 0.66 | 38 | 6 |
| training | Hallmark and GO Slim | 90 | 0.28 | 11 | 3 |
| acute | all collections | 1,370 | 0.35 | 287 | 53 |
| acute | collapse survivors | 95 | 0.58 | 95 | 19 |
| acute | Hallmark and GO Slim | 90 | 0.49 | 43 | 10 |

A set is discordant when its two NES differ in sign. Hallmark and GO Slim are paired because their
sets do not nest.

## Set classification clears chance only on the acute bout

singscore gives one rank-based score per set per sample on the imputed matrix. It carries no
p-value and never sees a contrast. Subject dominates raw scores: PC1 carries 25.6% of the variance,
of which subject explains 0.68, so step 05 reads within-subject change as well as levels.

Step 05 puts each score through the eight classification tasks of `02_Differential` and the
Spearman association in three windows, with BH within collection and task (classification) or
collection, window and outcome (association). Nominal hits over chance for classification:

| Task | Hallmark | KEGG_Legacy | Reactome | GOBP | GO_Slim |
|---|---:|---:|---:|---:|---:|
| training, HR | 0.59 | 1.29 | 0.79 | 0.60 | 1.72 |
| training, LR | 0.00 | 0.65 | 0.12 | 0.56 | 0.34 |
| acute, HR | 2.35 | 1.94 | 1.34 | 1.71 | 1.38 |
| acute, LR | 2.94 | 4.19 | 1.77 | 2.23 | 2.07 |
| HR vs LR at T1 *(floor)* | 0.59 | 0.00 | 0.12 | 0.24 | 0.00 |
| HR vs LR at T2 | 0.00 | 0.00 | 0.06 | 0.29 | 0.00 |
| HR vs LR, training change | 1.18 | 1.94 | 0.30 | 0.47 | 1.03 |
| HR vs LR, acute change | 1.18 | 0.97 | 0.37 | 0.98 | 0.34 |

The acute bout clears chance in every collection in both arms; no other task does so across
collections. No set survives BH in any task. The paired Wilcoxon has a floor set by the number of
pairs: 6 pairs (HR training) cannot go below p = 0.031, 7 (HR acute) below 0.016, 8 (LR) below
0.0078.

Four set-outcome pairs survive BH within collection: peroxisomal protein import (Reactome) and
peroxisome (KEGG) against `volume_load` over training, GO Slim DNA replication against baseline
`d_1rm_ext`, and Reactome basal-body anchoring against `comp_hypertrophy` over the acute bout.
Baseline level against `d_1rm_ext` runs at 3.4 times chance across collections.

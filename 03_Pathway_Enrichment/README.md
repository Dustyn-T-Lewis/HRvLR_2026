# 03 · Pathway Enrichment

The pathway level. Tests every gene set on the nine contrasts, scores every sample on every set,
and puts the scores through the protein level's classification and association.

| Step | Runs | Writes |
|---|---|---|
| [`00_Gene_Sets`](00_Gene_Sets/README.md) | frozen MSigDB snapshot, protein-to-gene map, size filter, GO Slim sets | `gene_sets.rds` |
| [`01_Scores`](01_Scores/README.md) | singscore per set and sample | `singscore.rds` |
| [`02_Contrasts`](02_Contrasts/README.md) | fgsea and fry per contrast, `collapsePathways` | `set_tests.rds` |
| [`03_Classify`](03_Classify/README.md) | AUC and Wilcoxon p on eight tasks, an ROC curve per nominal set | workbook and PDF |
| [`04_Associate`](04_Associate/README.md) | sample-level limma model per trait, change-score Spearman | workbook and PDF |
| [`05_Volcano`](05_Volcano/README.md) | protein volcanoes with pathway rings | 8 volcanoes |
| [`06_NES_Scatter`](06_NES_Scatter/README.md) | HR NES against LR NES, training and acute | 4 composites |

Every collection is its own screen: BH runs within collection.

## 1,378 sets pass the size filter

The first run fetches four MSigDB collections (release 2026.1.Hs) and writes an RDS with an md5 to
`00_Gene_Sets/c_data/cache/`. Later runs verify the checksum and need no network. Restore the
RDS and its `.md5` together; move both aside to rebuild. `goslim_generic.obo` sits beside it with
its own `.md5`; the script does not fetch it. `00_Gene_Sets/README.md` has its source and refresh
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

## No set classifier passes BH

singscore gives one rank-based score per set per sample on the imputed matrix. It carries no
p-value and never sees a contrast. Subject dominates raw scores: PC1 carries 25.6% of the variance,
of which subject explains 0.68.

Sets at p < 0.05 over sets tested, per collection. None reaches BH < 0.05 in any task.

| Task | Hallmark | KEGG_Legacy | Reactome | GOBP | GO_Slim |
|---|---:|---:|---:|---:|---:|
| Training, HR (T1 to T2) | 1 / 33 | 4 / 62 | 13 / 328 | 27 / 897 | 5 / 58 |
| Training, LR (T1 to T2) | 0 / 33 | 2 / 62 | 2 / 328 | 25 / 897 | 1 / 58 |
| Acute bout, HR (T2 to T3) | 4 / 33 | 6 / 62 | 22 / 328 | 77 / 897 | 4 / 58 |
| Acute bout, LR (T2 to T3) | 5 / 33 | 13 / 62 | 29 / 328 | 100 / 897 | 6 / 58 |
| HR vs LR at T1 | 1 / 33 | 0 / 62 | 2 / 328 | 11 / 897 | 0 / 58 |
| HR vs LR at T2 | 0 / 33 | 0 / 62 | 1 / 328 | 13 / 897 | 0 / 58 |
| HR vs LR, training change | 2 / 33 | 6 / 62 | 5 / 328 | 21 / 897 | 3 / 58 |
| HR vs LR, acute change | 2 / 33 | 3 / 62 | 6 / 328 | 44 / 897 | 1 / 58 |

HR's smallest attainable p is 0.031 (training, 6 pairs) and 0.016 (acute, 7); LR's is 0.0078 (8).

## The acute response tracks the training shift in fibre type

No term of the sample model has a BH hit. Between-person mCSA has the most sets at p < 0.05, 140
of 1,378, as at the protein level; between-person type I fibre area follows with 114.

Change-score pairs at BH < 0.05, BH within window, outcome and collection:

| Window | Outcome | Collection | Sets |
|---|---|---|---:|
| acute | pct_type1_share | Reactome | 85 |
| acute | pct_type1_share | GOBP | 3 |
| acute | d_type1_share | GOBP | 1 |
| acute | comp_hypertrophy | Reactome | 1 (basal-body anchoring) |
| acute | pct_leg_ext_1rm | GO_Slim | 1 (tRNA metabolism) |
| baseline | d_leg_ext_1rm | GO_Slim | 1 (DNA replication) |
| baseline | pct_leg_press_1rm | GO_Slim | 1 (ECM organisation) |
| training | volume_load_total_kg | Reactome, KEGG_Legacy | 2 (peroxisomal import, peroxisome) |

The 85 Reactome sets are mostly proteasome, ubiquitin, cell-cycle and NF-κB signalling. For all
but three, the acute-bout score change rises with the percent change in type I fibre share over
training (15 subjects, |rho| 0.63 to 0.83). The share comes from MyoVision counts of 142 to 1,119
fibres per biopsy. One subject (HR_S29) moved +105%; Spearman reads ranks, so that value weighs as
one rank. Each of the 20 outcomes, 3 windows and 5 collections carries its own BH, so these hits
are not corrected across outcomes.

# 03 · Pathway Enrichment

The pathway level. Tests every gene set on the nine contrasts, scores every sample on every set, and
puts the scores through the protein level's classification and association.

| Step | Runs | Writes |
|---|---|---|
| [`00_Gene_Sets`](00_Gene_Sets/README.md) | frozen MSigDB snapshot, GO Slim sets, protein-to-gene map, size filter | `gene_sets.rds` |
| [`01_Scores`](01_Scores/README.md) | singscore per set and sample | `singscore.rds` |
| [`02_Contrasts`](02_Contrasts/README.md) | fgsea and fry per contrast, `collapsePathways()`, ANCOVA on singscore | `set_tests.rds` |
| [`03_Classify`](03_Classify/README.md) | Wilcoxon p, AUC and median paired change on eight tasks | workbook, PDF |
| [`04_Associate`](04_Associate/README.md) | sample-level limma model per trait, change-score Spearman, arm slope difference, two-predictor model | workbook, PDF |
| [`05_Volcano`](05_Volcano/README.md) | protein volcanoes with pathway rings | workbook, PDF |
| [`06_NES_Scatter`](06_NES_Scatter/README.md) | HR NES against LR NES, training and acute | workbook, PDF |

Every collection (Hallmark, KEGG_Legacy, Reactome, GOBP, GO_Slim) is its own screen: BH runs within
collection for classification and association. Set contrasts pool the five collections within each
contrast.

Results live in each step workbook and PDF; each sub-stage README lists the sheets.
`05_Summary/02_Candidates` gathers the nominal hits across levels.

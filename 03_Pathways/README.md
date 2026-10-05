# 03 · Pathways

The pathway level. Tests every gene set on the nine contrasts, scores every sample on every set, and
puts the scores through the protein level's classification and association.

| Step | Runs | Writes |
|---|---|---|
| [`01_Gene_Sets`](01_Gene_Sets/README.md) | frozen MSigDB snapshot, GO Slim sets, protein-to-gene map, size filter | `gene_sets.rds` |
| [`02_Scores`](02_Scores/README.md) | singscore per set and sample | `singscore.rds` |
| [`03_Contrasts`](03_Contrasts/README.md) | fgsea, cameraPR and fry per contrast, `collapsePathways()`, ANCOVA on singscore | `set_tests.rds` |
| [`04_Classify`](04_Classify/README.md) | Wilcoxon p, AUC and median paired change on eight tasks | workbook, PDF |
| [`05_Associate`](05_Associate/README.md) | sample model, continuous change scores with their chance line, within-arm, arm slope, two-predictor, interaction, flips | workbook, PDF |

The volcano and NES-scatter pages moved to `06_Figures`.

Every collection (Hallmark, KEGG_Legacy, Reactome, GOBP, GO_Slim) is its own screen: BH runs within
collection for classification and association. Set contrasts pool the five collections within each
contrast.

Results live in each step workbook and PDF; each sub-stage README lists the sheets.
`05_Integration/02_Candidates` gathers the nominal hits across levels.

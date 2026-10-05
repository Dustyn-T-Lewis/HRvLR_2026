# 04 · Network

The module level. Finds proteins that move together, runs the protein level's contrasts,
classification and association on the module eigengenes, then names the modules and asks whether the
two arms share them.

| Step | Runs | Writes |
|---|---|---|
| [`01_Build`](01_Build/README.md) | WGCNA on subject-centred abundance | `modules.rds` |
| [`02_Contrasts`](02_Contrasts/README.md) | nine contrasts on eigengenes and with fry, ANCOVA, membership against protein results | workbook, PDF |
| [`03_Classify`](03_Classify/README.md) | Wilcoxon p, AUC and median paired change on eight tasks | workbook, PDF |
| [`04_Associate`](04_Associate/README.md) | sample-level limma model per trait, change-score Spearman, arm slope difference, two-predictor model | workbook, PDF |
| [`05_Characterise`](05_Characterise/README.md) | ORA against the stage 03 sets, hubs, STRING enrichment | workbook, PDF |
| [`06_Preserve`](06_Preserve/README.md) | modules built in each arm, preservation in the other | workbook, PDF |

Modules are defined on abundance centred within subject and scored on raw abundance, so they
describe co-movement rather than who a biopsy came from. Eigengene tests use the protein design and
subject block.

Results live in each step workbook and PDF; each sub-stage README lists the sheets.
`05_Integration/02_Candidates` gathers the nominal hits across levels.

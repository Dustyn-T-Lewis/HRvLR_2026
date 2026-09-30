# 02 · Differential Expression

The protein level. Tests every protein on nine contrasts, scores how well it separates the groups on
eight tasks, and ties it to phenotype per biopsy and per subject.

| Step | Runs | Writes |
|---|---|---|
| [`01_Design`](01_Design/README.md) | design formula, nine contrasts, within-subject correlation | `design.rds` |
| [`02_Contrasts`](02_Contrasts/README.md) | proteoDA limma fit, ANCOVA, `treat()`, fry concordance, unadjusted refit | `fit.rds` |
| [`03_Classify`](03_Classify/README.md) | Wilcoxon p, AUC and median paired change on eight tasks | workbook, PDF |
| [`04_Associate`](04_Associate/README.md) | sample-level limma model per trait, change-score Spearman, arm slope difference, two-predictor model | `associate.rds` |

Every protein test reads the unimputed normalised matrix, so limma fits each protein on the samples
where it was seen. BH runs within each contrast, task, trait and term, or window and outcome. The
pathway and module levels run the same three steps in the same slots.

Results live in each step workbook and PDF; each sub-stage README lists the sheets.
`05_Summary/01_Candidates` gathers the nominal hits across levels.

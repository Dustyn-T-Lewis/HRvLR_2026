# 02 · Differential Expression

The protein level. Tests every protein on nine contrasts, scores how well it separates the groups on
eight tasks, and ties it to phenotype per biopsy and per subject.

| Step | Runs | Writes |
|---|---|---|
| [`01_Design`](01_Design/README.md) | design formula, nine contrasts, within-subject correlation | `design.rds` |
| [`02_Contrasts`](02_Contrasts/README.md) | limpa fit (`dpcDE`), moderated ANCOVA, `treat()`, fry concordance, fit without blood, leave one subject out | `fit.rds` |
| [`03_Classify`](03_Classify/README.md) | Wilcoxon p, AUC and median paired change on eight tasks | workbook, PDF |
| [`04_Associate`](04_Associate/README.md) | sample-level limma model per trait, change-score Spearman, arm slope difference, two-predictor model | `associate.rds` |

The contrasts read limpa's quantification with its standard errors; classification and association
read the blood-corrected values where observed. BH runs within each contrast, task, trait and term, or window and outcome. The
pathway and module levels run the same three steps in the same slots.

Results live in each step workbook and PDF; each sub-stage README lists the sheets.
`05_Summary/02_Candidates` gathers the nominal hits across levels.

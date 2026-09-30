# 04_Associate

Ties each module eigengene to phenotype, the protein level's two ways.

| | |
|---|---|
| Reads | `01_Build_Modules/c_data/modules.rds`, `00_Input/phenotype.csv` |
| Writes | `c_data/04_associate.xlsx`, `b_reports/04_associate_figures.pdf` |
| Run | `Rscript 04_Network/04_Associate/a_script/04_associate.R` |

## Method

Sample model. Each T1 and T2 biopsy carries its own values for six traits: fibre area (mixed, type
I, type II), mCSA and the two 1RMs. Per trait, limma fits
`eigengene ~ timepoint + between + within` with subject blocked, then `eBayes(robust = TRUE)`:
`between` is the subject's mean trait, `within` the biopsy's deviation from it. Arm stays out, since
it was defined from fibre-area change. BH within trait and term.

Change score. One value per subject in three windows (training change, T1 level, acute change),
Spearman against 8 outcomes: `comp_hypertrophy`, `volume_load_total_kg`, and the percent change of
each trait. Pooled p uses the t approximation of `cor.test(exact = FALSE)`; within-arm p
is exact, from every ordering of the ranks. BH within window and outcome. Arm is the median split of
`comp_hypertrophy`, so pooled correlations with it and with fibre-area change partly carry the group
split.

Two limma models read the same per-subject values. Arm difference: `eigengene ~ arm * outcome`, whose
interaction is HR's slope minus LR's per outcome unit, the test of whether a eigengene tracks an
outcome differently in the two arms. Two predictor: `eigengene ~ comp_hypertrophy + strength gain`,
once per 1RM test, so growth and strength each come adjusted for the other (the two are unrelated
here). The two-thirds floor applies. BH within window and outcome, or window, strength test and
term.

## Figures

Overview pages: within- and between-person t for every trait beside status maps, then pooled rho for
every outcome, one page per window. Then every model hit at p < 0.05 (trait at the biopsy against
the module's value, one line per subject from T1 to T2) every change-score pair at p < 0.05 (one
point per subject, a least-squares line per arm drawn as a guide), twelve panels to a page, ordered
by p, then |t| or |rho|, then ID. Then every
arm-difference and two-predictor hit at p < 0.05, drawn the same way as the change scores.

## Workbook

| Sheet | Holds |
|---|---|
| `counts` | Tested, p < 0.05 and BH < 0.05 counts per term and trait, and per window and outcome. |
| `model` | Sample model: effect per trait unit, t, p, BH FDR and Storey q per module, trait and term. |
| `change_score` | Change score: pooled Spearman rho, p, BH FDR and Storey q per module, window and outcome. |
| `change_score_by_arm` | Change score within each arm, exact p, every module. |
| `arm_difference` | HR slope minus LR slope, t, p, BH FDR and Storey q per module, window and outcome. |
| `two_predictor` | Growth and strength terms, t, p, BH FDR and Storey q per module, window and test. |
| `plot_index` | PDF page of each overview, and page and panel position of every p < 0.05 row. |

# 04_Associate

Ties each protein to phenotype two ways: a sample-level model per biopsy and a change-score
correlation per subject.

| | |
|---|---|
| Reads | `01_Preprocess/02_Normalization/c_data/DAList_normalized.rds`, `01_Preprocess/03_Imputation/c_data/DAList_imputed.rds` (overview row order), `00_Input/phenotype.csv` |
| Writes | `c_data/associate.rds`, `c_data/04_associate.xlsx`, `b_reports/04_associate_figures.pdf` |
| Run | `Rscript 02_Differential_Expression/04_Associate/a_script/04_associate.R` |

## Method

Sample model. Each T1 and T2 biopsy carries its own values for nine traits: fibre area (mixed, type
I, type II), fibre count (mixed, type I), type I share, mCSA and the two 1RMs. Per trait, limma fits
`abundance ~ timepoint + between + within` with subject blocked, then `eBayes(robust = TRUE)`:
`between` is the subject's mean trait, `within` the biopsy's deviation from it. Arm stays out, since
it was defined from fibre-area change. A protein enters a trait's model, or a window's correlation,
only when two thirds of the biopsies or subjects observed it. BH within trait and term.

Change score. One value per subject in three windows (training change, T1 level, acute change),
Spearman against 20 outcomes: `comp_hypertrophy`, `volume_load_total_kg`, and the change and percent
change of each trait. Pooled p uses the t approximation of `cor.test(exact = FALSE)`; within-arm p
is exact, from every ordering of the ranks. BH within window and outcome. Arm is the median split of
`comp_hypertrophy`, so pooled correlations with it and with fibre-area change partly carry the group
split. `associate.rds` holds both tables; `04_Network/02_Contrasts` reads the training change
scores.

## Figures

Overview pages: within- and between-person t for every trait beside status maps, then pooled rho for
every outcome, one page per window. Then every model hit at p < 0.05 (trait at the biopsy against
the protein's value, one line per subject from T1 to T2) and every change-score pair at p < 0.05
(one point per subject, a least-squares line per arm drawn as a guide), twelve panels to a page,
ordered by p, then |t| or |rho|, then ID.

## Workbook

| Sheet | Holds |
|---|---|
| `counts` | Tested, p < 0.05 and BH < 0.05 counts per term and trait, and per window and outcome. |
| `model` | Sample model: effect per trait unit, t, p and BH FDR per protein, trait and term. |
| `change_score` | Change score: pooled Spearman rho, p and BH FDR per protein, window and outcome. |
| `change_score_by_arm` | Change score within each arm, exact p, rows with p < 0.05. |
| `plot_index` | PDF page of each overview, and page and panel position of every p < 0.05 row. |

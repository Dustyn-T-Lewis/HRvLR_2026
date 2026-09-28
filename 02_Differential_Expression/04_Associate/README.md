# 04_Associate

Ties each protein to phenotype two ways: a sample-level model per biopsy and a change-score
correlation per subject.

| | |
|---|---|
| Reads | `01_Preprocess/02_Normalization/c_data/DAList_normalized.rds`, `01_Preprocess/03_Imputation/c_data/DAList_imputed.rds` (overview row order), `00_Input/phenotype.csv` |
| Writes | `c_data/associate.rds`, `c_data/04_associate.xlsx`, `b_reports/04_associate_figures.pdf` |
| Run | `Rscript 02_Differential_Expression/04_Associate/a_script/04_associate.R` |
| Cost | about 3 min |

## Sample model

Each T1 and T2 biopsy carries its own values for nine traits: fibre area (mixed, type I, type II),
fibre counts (mixed, type I), type I share (type I count over mixed count), mCSA and the two 1RMs.
Per trait, limma fits `abundance ~ timepoint + between + within` with subject blocked by
`duplicateCorrelation()`, then `eBayes(robust = TRUE)`. `between` is the subject's mean trait and
`within` the biopsy's deviation from it. Arm stays out, since it was defined from fibre-area
change. A protein enters only when two thirds of the trait's biopsies saw it. Effects are per trait
unit; BH runs within trait and term.

## Change score

One value per subject in three windows: training change (T2 − T1), T1 level and acute change
(T3 − T2). Spearman against 20 outcomes: `comp_hypertrophy`, `volume_load_total_kg`, and the change
and percent change of each of the nine traits. A protein is tested only when two thirds of the
window's subjects saw it. Pooled p is the t approximation `cor.test(exact = FALSE)` uses; above
nine observations `cor.test`'s default p is an Edgeworth series that returns 0 in the tail. Within
each arm (4 to 8 subjects) p is exact, from every ordering of the ranks. BH within window and
outcome.

## Outputs

`associate.rds` holds `change_score` and `model`; `04_Network/02_Contrasts` reads the training
change scores. `counts` gives tested, p < 0.05 and BH < 0.05 counts for every model term and trait
and every window and outcome.

The figure PDF opens with overview pages: one for the sample model (within- and between-person t for
every trait, beside their status maps) and one per change-score window (pooled rho for every
outcome, beside its status map). Rows follow average-linkage clustering (1 - Pearson r) of the
imputed matrix, for display only. The grids follow. Every model hit at p < 0.05 gets a panel: trait
at the biopsy against abundance, one line per subject from T1 to T2, and a header with the effect,
p, BH q and the biopsies and subjects where the protein was observed. Every change-score pair at p <
0.05 gets a scatter, one point per subject, with a least-squares line per arm drawn as a guide; the
tests are Spearman. The header gives the pooled rho, p, BH q and subjects used, the corner each
arm's rho and exact p. Twelve panels to a page, ordered by p, then |t| or |rho|, then ID.
`plot_index` gives the PDF page and panel of every drawn row. Arm is the median split of
`comp_hypertrophy`, so pooled correlations with it and with fibre-area change partly carry the group
split; the change-score captions say so.

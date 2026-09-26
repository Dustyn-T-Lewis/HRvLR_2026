# 04_Associate

Ties each protein to phenotype two ways: a sample-level model per biopsy and a change-score
correlation per subject.

| | |
|---|---|
| Reads | `01_Preprocess/02_Normalization/c_data/DAList_normalized.rds`, `00_Input/phenotype.csv` |
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
change scores. `chance_expectation` gives nominal counts against 5% of tested for every model term
and trait and every window and outcome.

The figure PDF opens on the chance ratios, then gives every nominal model hit a panel (one line per
subject from T1 to T2, coloured by arm), then every nominal change-score pair a scatter
(`ggpubr::ggscatter`, a line per arm, each arm's rho and exact p in the corner). Twelve panels to a
page, by p. Each caption gives the tested count and the chance expectation.

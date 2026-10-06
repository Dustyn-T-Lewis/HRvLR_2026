# F3_People

Figure 3: can one feature tell HR from LR, and does a protein's link to strength gain differ
between them?

| | |
|---|---|
| Reads | `02_HRvLR_groups/05_Classify/{01_Proteins,02_Sets,03_Modules}/c_data/*.xlsx` (`auc`, `counts`), `02_HRvLR_groups/06_Associate/c_data/06_associate.xlsx` (`group_difference`, `within_group`), `04_Continuous_no_group/01_Associate/01_Proteins/c_data/01_proteins.xlsx` (`interaction`), `01_Trunk/03_Quantify/c_data/DAList_corrected.rds`, `00_Input/phenotype.csv` |
| Writes | `b_reports/f3_people.pdf` (300 x 175 mm), `c_data/f3_people.xlsx` |
| Run | `Rscript 02_HRvLR_groups/Figures/F3_People/a_script/f3_people.R` |

The script reads saved outputs only. Its one computation is the Spearman correlation across
proteins in panel b, a description with no test. Before writing, it checks the classification
counts against their sources and recomputes panel d's within-group rho from the plotted points.

## Panels

| Panel | Shows | Rule |
|---|---|---|
| a | The best single protein, set and module for each of the five between-group tasks: AUC with its 95% interval (pROC's stratified bootstrap where the stage drew the feature, DeLong's otherwise), and beside it the features at p < 0.05 against the number chance gives. | smallest Wilcoxon p per task and level |
| b | Every protein's group-by-trait difference (t, HR slope minus LR slope) against its continuous growth-by-strength interaction (t), for the window and trait of the lead group difference. | lead: smallest group-difference p among proteins |
| c | The same proteins' Spearman rho with the trait in HR against LR; group-difference hits at FDR < 0.05 marked. | as b |
| d | The lead protein's change against the trait, person by person, with a line and Spearman rho per group. | as b |

Single features do not separate the groups beyond chance. The best protein per task reaches AUC
0.93 to 0.98, but the number of proteins at p < 0.05 sits at or below the 95 that chance gives
(44 to 100), and no module reaches p < 0.05 on any task. A best-of-1,900 AUC near 1 with seven
and eight people is what a null screen produces.

The group differences that reach FDR < 0.05 all sit in one window and trait: training change
against leg press 1RM % change (11 proteins, 28 sets, 1 module in `06_Associate`). HR and LR are a
median split of the hypertrophy composite, so the group-by-trait test is a dichotomised version of
the continuous growth-by-strength interaction; across all 1,900 proteins the two t statistics have
Spearman 0.88, and 11 of the 20 continuous hits at FDR < 0.05 are the same proteins. The hits are
sarcomere proteins (TPM3, MYH1, MYH7, MYL2, MYBPC2, TTN, TNNT1 and others); in every one the
within-group rho is negative in HR and positive in LR. Within-group rho rests on six to eight
people and is descriptive; the group difference is the test.

Module preservation between the groups is in `04_Modules/02_Preserve` and not drawn here.

## Workbook

| Sheet | Holds |
|---|---|
| `panel_sources` | Per panel: what it shows, the file and sheet it reads, and the rule that picks its rows. |
| `panel_a` | Best feature per task and level: AUC, DeLong and bootstrap intervals, p, counts at p < 0.05. |
| `panel_b` | Per protein: group-difference t and continuous growth x strength t, with p and FDR. |
| `panel_c` | Per protein: Spearman rho in each group, and the group-difference evidence. |
| `panel_d` | The lead protein's training change and the trait, per person. |

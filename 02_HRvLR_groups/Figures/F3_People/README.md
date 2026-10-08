# F3_People

Parked draft: the script runs and its outputs are rebuilt, but the figure is not final.

Figure 3: can one feature tell the groups apart, and does a protein's link to strength gain differ between them? Panel a reads `06_Classification`; panels b to d follow the protein with the smallest group-by-trait p in `07_Trait_Associations`.

| Panel | Draws | Rule that picks the rows |
|---|---|---|
| a | For each level (proteins, sets, modules) and each of the five classification tasks, the best single feature's AUC with its 95% interval: pROC's stratified bootstrap where the classification stage drew the feature, DeLong's otherwise (solid and dotted lines). Beside it, the feature's name and the count of features at p < 0.05 against the count chance gives (5% of those tested). | smallest Wilcoxon p per task and level, then AUC furthest from 0.5 |
| b | Every protein's group-difference t (HR slope minus LR slope, from `feature ~ group * trait`) against its t for the continuous growth-by-strength interaction from the pooled analysis. The subtitle gives the Spearman correlation of the two across proteins, a description. Up to eight group-difference hits at FDR < 0.05 are named. | the lead protein's window and trait |
| c | The same proteins' within-group Spearman rho with the trait, HR on x and LR on y, with group-difference hits at FDR < 0.05 marked and named; the subtitle counts how many of them change sign between groups. | as b |
| d | The lead protein's training change against the trait, one point per person, a least-squares line per group, each group's rho in the legend and the group-difference t, p and FDR in the subtitle. | the lead: smallest group-difference p among proteins |

| | |
|---|---|
| Reads | `02_HRvLR_groups/06_Classification/01_Proteins/c_data/protein_classification.xlsx`, `02_HRvLR_groups/06_Classification/02_Sets/c_data/set_classification.xlsx`, `02_HRvLR_groups/06_Classification/03_Modules/c_data/module_classification.xlsx` (each `auc`, `counts`), `02_HRvLR_groups/07_Trait_Associations/c_data/trait_associations.xlsx` (`group_difference`, `within_group`), `04_Continuous_no_group/01_Trait_Associations/01_Proteins/c_data/protein_trait_associations.xlsx` (`interaction`), `01_Processing/03_Quantification/c_data/DAList_corrected.rds`, `00_Input/phenotype.csv` |
| Writes | `b_reports/f3_people.pdf`, `c_data/f3_people.xlsx` |
| Run | `Rscript 02_HRvLR_groups/Figures/F3_People/a_script/f3_people.R` |

## Read with care

- The script completes only when the lead protein's window is `training` and its trait a 1RM % change: the continuous interaction in panel b covers only those traits, and panel d draws the T2 - T1 change. Otherwise it stops rather than draw a wrong figure.
- Panel a's best feature per task is the extreme of hundreds of tests; the count against chance beside it is the fairer summary.
- The lead protein is chosen for the smallest p of many, so its group difference is inflated; panels b and c place it among all proteins.
- Within-group rho rests on six to eight people and is descriptive.

## Sheets

| Sheet | One row is | Key columns |
|---|---|---|
| `panel_sources` | a panel | `panel`, `shows`, `file`, `sheet`, `rule` |
| `panel_a` | a level and task | `feature`, `label`, `auc`, `auc_lower`, `auc_upper`, `boot_lower`, `boot_upper`, `p`, `fdr`, `n_positive`, `n_negative`, `n_nominal`, `n_tested` |
| `panel_b` | a protein | `t_group`, `p_group`, `fdr_group`, `t_continuous`, `p_continuous`, `fdr_continuous`, `hit` |
| `panel_c` | a protein | `rho_HR`, `rho_LR`, `hit` |
| `panel_d` | a person | `subject`, `value`, `trait`, `arm` |

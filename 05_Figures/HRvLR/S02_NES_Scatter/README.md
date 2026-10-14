# S02_NES_Scatter

Supplement: do HR and LR shift the same gene sets, set by set? Each set's fgsea NES in one group is plotted against its NES in the other, for training (`Training_HR` against `Training_LR`) and for the bout (`Acute_HR` against `Acute_LR`). A set on the dashed identity line moved the same way and as far in both groups.

| Term | Meaning |
|---|---|
| significance | where the set reaches fgsea FDR < 0.05: both, one contrast, or neither (NS) |
| discordant | the two NES have opposite signs; the discordant panel and counts take only sets significant in at least one group |
| survivor | kept as non-redundant by `collapsePathways()` in either contrast |
| concordance rho | Spearman correlation of the two NES across the sets in a population |

| | |
|---|---|
| Reads | `03_Group_comparisons/HRvLR/c_data/03_set_tests.rds` (`set_tests`) |
| Writes | `b_reports/s02_nes_scatter.pdf`, `c_data/s02_nes_scatter.xlsx` |
| Run | `Rscript 05_Figures/HRvLR/S02_NES_Scatter/a_script/s02_nes_scatter.R` |

## Read with care

- Sets overlap and share proteins; the concordance rho describes agreement, and its p overstates the evidence. The Hallmark and GO Slim panels limit the nesting.
- NES compares a set with the other proteins in the same contrast, so a group whose whole proteome moved more can differ in NES without the set behaving differently.
- Agreement between the groups' NES is not a test of a group difference; the interaction contrasts in `03_set_contrasts` are.

## Sheets

| Sheet | One row is | Key columns |
|---|---|---|
| `concordance` | a pair and population | `pair`, `population`, `n_sets`, `n_significant`, `n_discordant`, `rho`, `p`, `same_sign` |
| `nes_scatter` | a set scored in both contrasts of a pair | `pair`, `x_contrast`, `y_contrast`, `collection`, `pathway`, `n_proteins`, `nes_x`, `nes_y`, `fdr_x`, `fdr_y`, `significance`, `discordant`, `survivor` |

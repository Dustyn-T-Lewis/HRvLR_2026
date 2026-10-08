# S03_Signatures

Supplement to Figure 1: does the signature carry-over hold with stricter sets, and where do the training-responsive proteins go after the bout? Panels a and b repeat F1's barcodes with sets chosen at Π < 0.05 (p^|logFC|). Panel c follows each group's p < 0.05 training proteins through T1, T2 and T3 in both groups, as the mean change from each person's T1.

| | Measured in the group that chose the set | Measured in the other group |
|---|---|---|
| Training window | the selection window: the set was chosen here, so it carries no test | fry on that group's training contrast |
| Bout window | shares T2 with the selection, so a protein high at T2 by chance reads as up then down; faded and marked † | fry on that group's acute contrast |

| | |
|---|---|
| Reads | `02_HRvLR_groups/02_Protein_Contrasts/c_data/protein_contrasts.xlsx` (`fry_concordance`, `set_trajectory`), `02_HRvLR_groups/03_Set_Contrasts/c_data/set_tests.rds` (`protein_results`) |
| Writes | `b_reports/s03_signatures.pdf`, `c_data/s03_signatures.xlsx` |
| Run | `Rscript 02_HRvLR_groups/Figures/S03_Signatures/a_script/s03_signatures.R` |

## Read with care

- Π sets are smaller than p < 0.05 sets, so a weaker fry result for a Π set does not contradict the p < 0.05 result.
- In panel c the own-group bout window shares T2 with the selection (marked †); the other-group cells are the clean tests.
- Panel c's standard error is across people and ignores that the choosing group's members were picked on the same data.

## Sheets

| Sheet | One row is | Key columns |
|---|---|---|
| `pi_barcode` | a lane and direction (panels a, b) | `window`, `chosen`, `ranked_in`, `direction`, `n`, `same_sign_share`, `fdr` |
| `trajectory` | a set, group and time point (panel c) | `set_from`, `set`, `n_proteins`, `arm`, `timepoint`, `mean_change`, `se` |
| `window_stats` | a set measured in a group over a window (panel c) | `chosen`, `arm`, `window`, `direction`, `n_proteins`, `same_sign_share`, `fdr` |

# S03_Signatures

Supplement to Figure 1: the signature ranks with the smaller Π sets, and the training sets followed
through the bout.

| | |
|---|---|
| Reads | `02_HRvLR_groups/02_Protein_Contrasts/c_data/protein_contrasts.xlsx` (`fry_concordance`, `set_trajectory`), `02_HRvLR_groups/03_Set_Contrasts/c_data/set_tests.rds` |
| Writes | `b_reports/s03_signatures.pdf` (300 x 110 mm), `c_data/s03_signatures.xlsx` |
| Run | `Rscript 02_HRvLR_groups/Figures/S03_Signatures/a_script/s03_signatures.R` |

## Panels

| Panel | Shows |
|---|---|
| a, b | F01 e and f with the Π < 0.05 sets (13, 7, 11 and 14 training proteins; 21, 12, 44 and 18 acute). |
| c | Each arm's p < 0.05 training sets followed through T1, T2 and T3 in both arms: mean change from each person's T1 (log2) ± SE, with the training window shaded blue and the bout orange. |

Panel c reports, in each window, the share of the set's proteins moving as they did in training and
fry's FDR on that window's contrast. Two cells need care. In the arm that chose a set, the training
window is where it was picked, so it carries no test ("selection window"). The bout window there
shares T2 with the selection: a protein high at T2 by chance reads as up, then down. Those lines are
faded and their bout result marked †. The other-arm cells are clean:

| Training sets | Measured in | Training window | Bout window |
|---|---|---|---|
| HR's, up / down | LR | 90% / 74% same way, FDR 0.14 / 0.25 | 27% / 46%, FDR 0.64 / 0.92 |
| LR's, up / down | HR | 82% / 76%, FDR 0.016 / 0.026 | 30% / 45%, FDR 0.13 / 0.57 |

Training-responsive proteins tend to drift back toward baseline after the bout, but no clean cell
reaches FDR 0.05 for that step.

## Workbook

| Sheet | Holds |
|---|---|
| `pi_barcode` | Panels a, b: per window and lane, set size, share moving the same way and fry FDR. |
| `trajectory` | Panel c: mean change from T1 and SE per set, arm and time point. |
| `window_stats` | Panel c: per cell and window, share moving as in training and fry FDR. |

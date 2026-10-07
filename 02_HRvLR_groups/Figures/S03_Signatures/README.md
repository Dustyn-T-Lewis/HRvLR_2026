# S03_Signatures

Supplement to Figure 1: does the signature carry-over hold with stricter sets, and where do the
training-responsive proteins go after the bout?

## The idea

F1 e and f take each group's proteins at p < 0.05 in a window and ask, with fry, whether they move
the same way in the other group. Panels a and b repeat that with the sets chosen at
Π < 0.05. Π is p raised to the power \|logFC\| (Xiao et al. 2014), so it keeps proteins with both
a small p and a large fold change. The barcode is drawn as in F1: members as ticks at their rank
in the other group's moderated t, up set above the axis, down set below, with a tricube moving
average of membership; each lane gives set size, the share moving the same way and fry's FDR.

Panel c follows each group's p < 0.05 training proteins through T1, T2 and T3, in both groups. For
each set (up and down, chosen in HR or in LR) and each biopsy, the members' mean log2 abundance is
taken; each person's T1 value is subtracted; the mean and standard error across people are drawn
per group and time point. The training window is shaded blue, the bout orange.

| | Measured in the group that chose the set | Measured in the other group |
|---|---|---|
| Training window | the selection window: the set was chosen here, so it carries no test | fry on that group's training contrast |
| Bout window | shares T2 with the selection, so a protein high at T2 by chance reads as up then down; faded and marked † | fry on that group's acute contrast |

Each window label gives the share of the set's proteins moving as they did in training and fry's
FDR on that window's contrast, from `fry_concordance` in `02_Protein_Contrasts`. The two other-group
cells are the clean ones: different people, no shared cell mean.

## What the code does

| Step | Function | Detail |
|---|---|---|
| Panels a, b | `signature_barcode(set_level = "pi")` in `R/figure_panels.R` | Rankings from `protein_results`; fry numbers from the `pi_` rows of `fry_concordance`. |
| Panel c | `ggplot2`, `ggh4x::facet_grid2()` | `set_trajectory` rows for the `nominal_` sets: a 2 x 2 grid, rows by the group that chose the set, columns by the group measured. Own-group lines at reduced opacity. |
| Window labels | `fry_concordance` rows from `Training_` sets | Other group's training window, and both groups' bout windows; the own group's training window reads "selection window". |
| Assembly | `patchwork`, `save_figure()` | 300 x 110 mm, `cairo_pdf`. |
| Checks | `stopifnot()` | Eight barcode lanes with fry numbers; 24 trajectory rows (four sets, two groups, three time points); twelve window results; two labels per cell and direction. |

## How to read the outputs

- Π < 0.05 is stricter than p < 0.05 for any protein with \|logFC\| below 1, so Π sets are usually
  smaller and fry has fewer proteins to work with. A weaker fry result for a Π set does not
  contradict the p < 0.05 result.
- The standard error in panel c is across people within a group and ignores that the members were
  chosen on the same data in the choosing group.
- People with no T1 biopsy drop out of panel c, since change is measured from each person's own T1.
- fry's FDR in each label runs over the sets tested together in one call in `02_Protein_Contrasts`,
  not over the whole figure.

## Files

| | |
|---|---|
| Reads | `02_HRvLR_groups/02_Protein_Contrasts/c_data/protein_contrasts.xlsx` (`fry_concordance`, `set_trajectory`), `02_HRvLR_groups/03_Set_Contrasts/c_data/set_tests.rds` (`protein_results`) |
| Writes | `b_reports/s03_signatures.pdf`, `c_data/s03_signatures.xlsx` |
| Run | `Rscript 02_HRvLR_groups/Figures/S03_Signatures/a_script/s03_signatures.R` |

```
S03_Signatures/
├── a_script/s03_signatures.R
├── b_reports/s03_signatures.pdf
└── c_data/s03_signatures.xlsx
```

## Workbook: `s03_signatures.xlsx`

`read_me` lists the sheets; `columns` defines every column; `input_manifest` and
`package_versions` close the workbook.

| Sheet | One row is | Key columns |
|---|---|---|
| `pi_barcode` | a lane and direction (panels a, b) | `window`, `chosen`, `ranked_in`, `direction`, `n`, `same_sign_share`, `fdr` |
| `trajectory` | a set, group and time point (panel c) | `set_from`, `set`, `n_proteins`, `arm`, `timepoint`, `mean_change`, `se` |
| `window_stats` | a set measured in a group over a window (panel c) | `chosen`, `arm`, `window`, `direction`, `n_proteins`, `same_sign_share`, `fdr` |

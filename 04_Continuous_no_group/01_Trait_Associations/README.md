# 01_Trait_Associations

Which proteins, gene-set scores and module eigengenes track each size and strength trait, with
every person pooled? The three levels call `associate_step()` in `R/associate.R`.

| Step | Values | Differs from the others |
|---|---|---|
| [`01_Proteins`](01_Proteins/README.md) | limpa's complete blood-corrected matrix | carries `observed_share`; saves `associate.rds` for `02_Trait_Pathways` and `03_Profile_Clusters` |
| [`02_Sets`](02_Sets/README.md) | singscore per qualifying gene set | BH within collection |
| [`03_Modules`](03_Modules/README.md) | module eigengenes | none |

BH runs within each screen:

| Analysis | One screen is |
|---|---|
| Sample model | one trait and term (`within`, `between`) |
| Change score | one window and trait |
| Two predictor | one window, 1RM and coefficient |
| Interaction | one window and 1RM |

The report draws up to 24 hits per screen; the workbook keeps every row, and `plot_index` gives
the page and panel of each drawn one.

# 05 · Figures

Five manuscript figures, each assembled from stage 01 to 04 outputs. F01 panel A also reads
`00_Input/phenotype.csv`. Diagnostics stay in each step's `b_reports/`.

| Figure | Shows | Panels |
|---|---|---|
| [`F01`](F01/README.md) | cohort and filtering | 3 |
| [`F02`](F02/README.md) | proteome | 4 |
| [`F03`](F03/README.md) | pathways | 3 |
| [`F04`](F04/README.md) | networks | 4 |
| [`F05`](F05/README.md) | classification and association across levels | 2 |

```sh
for f in F01 F02 F03 F04 F05; do Rscript 05_Figures/$f/a_script/$f.R; done
```

Each figure folder holds `a_script/panels/<letter>_<content>.R`, one script per panel, and
`a_script/<figure>.R`, which sources the panels and stitches them with patchwork at 178 mm wide.
Panels set titles, axes and colours; the composite sets layout, letters and size. `theme.R` holds
the shared theme, palettes and savers every script sources. The figures take under a minute.

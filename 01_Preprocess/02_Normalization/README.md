# 01_Preprocess / 02_Normalization

Compares normalisation methods and applies cyclic loess.

| | |
|---|---|
| Script | `a_script/02_normalize.qmd` |
| Reads | `01_Filtering/c_data/DAList_filtered.rds` |
| Writes | `c_data/DAList_normalized.rds`, `c_data/02_normalize.xlsx`, three proteoDA PDFs in `b_reports/` |

```sh
quarto render 01_Preprocess/02_Normalization/a_script/02_normalize.qmd --output-dir ../b_reports
```

`norm_comparison.pdf` draws each proteoDA method against Group_Time; `qc_pre.pdf` and
`qc_post.pdf` bracket the chosen one. The workbook holds the normalised matrix, the leading
components and each protein's eta squared on Group_Time.

## The span

proteoDA's cyclic loess is `limma::normalizeCyclicLoess(method = "fast")`. With
`adaptive.span = TRUE` the span comes from `chooseLowessSpan(nrow)`, not the 0.7 in the formals,
so it moves whenever filtering changes the protein count.

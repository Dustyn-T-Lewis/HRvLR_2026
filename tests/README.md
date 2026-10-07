# tests

`Rscript -e 'testthat::test_dir("tests/testthat")'` runs everything here in about 40 seconds;
`run_all.R` runs it after a full run.

| File | Checks |
|---|---|
| `test-storey_q.R` | Storey's q equals BH below 100 tests, qvalue's bootstrap estimate from 100 on, BH when qvalue fails |
| `test-classify_row.R` | AUC is W / (n1 n2), p is `wilcox.test`'s, the interval is pROC's DeLong; AUC below 0.5 kept; no interval at perfect separation |
| `test-associations.R` | `spearman_by_row()` against `cor.test` and `DescTools::CorCI`; `subject_fit()` slope against `lm` |
| `test-limpa.R` | `fit_limpa()` on limpa's simulated data; `blood_corrected()` against an `lm` slope |
| `test-grouping.R` | `relabel_samples()` and `for_design()` keep the grouped people and relabel them; per-person values and changes |
| `test-phenotype.R` | % change and the trait table each window meets |
| `test-validity.R` | the circular, caution and clean tags against the table in `03_Cluster_groups/Proteome/README.md`; `stack_inputs()` |
| `test-workbook.R`, `test-tables.R` | workbook sheet order, manifest md5, the column dictionary, page layout helpers |
| `test-outputs.R` | on the saved outputs: sheet order, every interval contains its estimate, contrast names follow the groups, the five classification tasks, validity on every stacked proteome row |

The unit tests cover every line of `R/helpers.R` and `03_Cluster_groups/Proteome/R/inputs.R`
(covr). The step functions (`classify_step()`, `associate_step()`, the contrast steps) are checked
through their outputs in `test-outputs.R`. Those checks skip when no outputs exist.

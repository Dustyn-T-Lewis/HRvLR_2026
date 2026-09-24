# HRvLR Proteomics

DIA mass-spectrometry proteomics of skeletal muscle from a resistance training study. Sixteen
subjects were biopsied at baseline (T1), after the training programme (T2) and after an acute bout
at the end of it (T3). They were labelled High or Low Responder by a median split of a composite
hypertrophy score, eight each. 48 MS runs; three are dropped as outliers, leaving 45.

The primary comparison is the training interaction: whether the proteome changed differently over
training in HR than in LR. The study also asks whether any protein, pathway or co-expression module
tracks how much a subject adapted, across ten phenotypes.

## Every stage passes data through disk

| Stage | Runs | Writes |
|---|---|---|
| [`00_Input/`](00_Input/README.md) | study data, two builders for derived inputs | `phenotype.csv`, `RBC_proteome_reference.tsv` |
| [`01_Preprocess/`](01_Preprocess/README.md) | filtering, cyclic loess, missForest | `DAList_normalized.rds`, `DAList_imputed.rds` |
| [`02_Differential_Expression/`](02_Differential_Expression/README.md) | design, nine contrasts, protein classification, protein against phenotype | `design.rds`, `fit.rds`, `phenotype.rds` |
| [`03_Pathway_Enrichment/`](03_Pathway_Enrichment/README.md) | gene set tests, per-sample set scores, set classification and association | `gene_sets.rds`, `set_tests.rds`, `singscore.rds` |
| [`04_Network/`](04_Network/README.md) | co-expression modules, preservation between arms, the same tests on modules | `modules.rds` |
| [`05_Figures/`](05_Figures/README.md) | five manuscript figures from stage 01 to 04 outputs | `F01.pdf` to `F05.pdf`, with PNGs and per-panel files |

Each sub-stage holds `a_script/` (code), `b_reports/` (`<step>_figures.pdf` where the step draws,
and the rendered HTML report for notebooks, which git ignores), `c_data/` (one workbook per step,
plus the `.rds` a later step reads) and a README. Every workbook opens on a
`read_me` sheet and ends with `input_manifest` and `package_versions`. Stages 01 and 02 are Quarto
notebooks; stages 03 and 04 are R scripts. Data passes through disk, so any sub-stage re-runs on its
own once its inputs exist.

## One command per step, in order

```sh
Rscript -e 'renv::restore()'

quarto render 01_Preprocess/01_Filtering/a_script/01_filter.qmd        --output-dir ../b_reports
quarto render 01_Preprocess/02_Normalization/a_script/02_normalize.qmd --output-dir ../b_reports
quarto render 01_Preprocess/03_Imputation/a_script/03_impute.qmd       --output-dir ../b_reports

quarto render 02_Differential_Expression/01_Design/a_script/01_design.qmd             --output-dir ../b_reports
quarto render 02_Differential_Expression/02_Differential/a_script/02_differential.qmd --output-dir ../b_reports
quarto render 02_Differential_Expression/03_Phenotype/a_script/03_phenotype.qmd       --output-dir ../b_reports

for s in 00_build_gene_sets 01_run_fgsea_and_fry 02_enrich_volcano_fgsea \
         03_enrich_scatter_fgsea 04_run_singscore 05_classify_and_associate_sets; do
  Rscript 03_Pathway_Enrichment/$s/a_script/$s.R
done

for s in 01_build_modules 02_characterise_modules 03_preserve_modules 04_test_modules \
         05_classify_and_associate_modules; do
  Rscript 04_Network/$s/a_script/$s.R
done

for f in F01 F02 F03 F04 F05; do Rscript 05_Figures/$f/a_script/$f.R; done
```

Stages 01 to 04 take about ten minutes. `04_Network/02` needs the STRING files;
`00_Input/README.md` has the download command.

## proteoDA fits the unimputed matrix; BH never pools across screens

Preprocessing and the fit use proteoDA: `DAList`, `zero_to_missing`, `filter_proteins_by_group`,
`filter_samples`, `normalize_data("cycloess")`, then `add_design`, `add_contrasts`,
`fit_limma_model` and `extract_DA_results`. The design is six cell means with subject as a random
effect.

The fitted matrix stays unimputed: limma fits each protein on the samples where it was seen. A
missForest copy is read only by the methods that need a complete matrix: `fry`, singscore and
WGCNA.

Every classification and association screen reports its nominal count beside the count chance
predicts. BH runs within each contrast, task, or window and outcome. Set tests pool the five
collections within a contrast; set classification and association split by collection.

## renv.lock pins every package

Stages 01 and 02: `proteoDA`, `limma`, `missForest`, `ranger`, `lme4`, `callr`,
`here`, `dplyr`, `tidyr`, `tibble`, `purrr`, `stringr`, `forcats`, `readr`, `readxl`, `ggplot2`,
`patchwork`, `writexl`, `sessioninfo`. Stage 03 adds `fgsea`, `singscore`, `msigdbr`, `GO.db`,
`GSEABase`, `AnnotationDbi`, `pROC`, `ggrepel` and `enrichVolcano`. Stage 04 adds `WGCNA`,
`clusterProfiler`, `enrichplot`, `STRINGdb`, `tidygraph`, `ggraph` and `igraph`. Rendering needs
Quarto.

`enrichVolcano` 0.3.0.9000 is recorded in `renv.lock` without a remote, so `renv::restore()`
cannot fetch it. To move to the GitHub release, install it and snapshot, then rerun stage 03 and
check its figures:

```r
renv::install("Dustyn-T-Lewis/enrichVolcano")
renv::snapshot()
```

The code calls only its exported `volcano_ring()`, `volcano_ring_theme()` and `ev_clean_label()`.

## Version 1 is the tag legacy-v1

The first version of this project, with its prediction screen (F06), `functions/`, `tests/` and
`archive/`, is the git tag `legacy-v1`.

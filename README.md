# HRvLR Proteomics

DIA mass-spectrometry proteomics of skeletal muscle from a resistance training study. Sixteen
subjects were biopsied at baseline (T1), after the training programme (T2) and after an acute bout
at the end of it (T3). They were labelled High or Low Responder by a median split of a composite
hypertrophy score, eight each. 48 MS runs; three are dropped as outliers, leaving 45.

The primary comparison is the training interaction: whether the proteome changed differently over
training in HR than in LR. The study also asks whether any protein, pathway or co-expression module
tracks how much a subject adapted, across ten phenotypes.

## Stages

| Stage | Contents | State |
|---|---|---|
| `00_Input/` | study data, two builders for derived inputs | ready |
| `01_Preprocess/` | protein report to a normalised DAList, and an imputed copy | ready |
| `02_Differential_Expression/` | model fitting, nine contrasts, protein classification, protein against phenotype | ready |
| `03_Pathway_Enrichment/` | gene set tests, per-sample set scores, set classification and phenotype association | ready |
| `04_Network/` | co-expression modules, what they are, and the same tests on them | ready |
| `05_Figures/` | manuscript panels | planned |

Each sub-stage holds `a_script/` (code), `b_reports/` (HTML reports or figures) and `c_data/`
(outputs), with a README. Stages 01 and 02 are Quarto notebooks; stages 03 and 04 are plain R
scripts. Data passes through disk, so any sub-stage re-runs on its own. "Planned" means only a
README exists.

## Data

The protein report, sample sheet and phenotype table are committed in `00_Input/`. The STRING v12
files `04_Network` reads are too large for git; `00_Input/README.md` has the download command.

## Running the pipeline

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

for s in 01_build_modules 02_characterise_modules 03_classify_and_associate_modules; do
  Rscript 04_Network/$s/a_script/$s.R
done
```

Everything above takes about five minutes.

## Approach

Preprocessing and the fit use **proteoDA**: `DAList`, `zero_to_missing`, `filter_proteins_by_group`,
`filter_samples`, `normalize_data("cycloess")`, then `add_design`, `add_contrasts`,
`fit_limma_model` and `extract_DA_results`. The design is six cell means with subject as a random
effect, because HR and LR are different people and a fixed subject term would absorb every
between-arm contrast.

The fitted matrix stays unimputed: limma fits each protein on the samples where it was seen. A
missForest copy is read only by the three methods that need a complete matrix: `fry`, singscore and
WGCNA.

Every screen reports its nominal count beside the count chance predicts, and BH runs within each
contrast, task, window or collection, never across them.

## Dependencies

Pinned in `renv.lock`. `proteoDA`, `limma`, `missForest`, `here`, `dplyr`, `tidyr`, `tibble`,
`purrr`, `ggplot2`, `patchwork`, `writexl`, `readxl`. Stage 03 adds `fgsea`, `singscore`, `msigdbr`,
`GO.db`, `GSEABase`, `AnnotationDbi`, `pROC`, `ggrepel`, `qpdf` and `enrichVolcano` (not on CRAN;
`renv::hydrate("enrichVolcano")` links a local install). Stage 04 adds `WGCNA`, `lme4`,
`clusterProfiler`, `enrichplot` and `STRINGdb`. Rendering needs Quarto.

Earlier designs of this project, a continuous phenotype sweep and a blind subtype search, are in
git history.

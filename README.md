# HRvLR Proteomics

DIA mass-spectrometry proteomics of skeletal muscle from a resistance training study. Sixteen
subjects were biopsied at baseline (T1), after the training programme (T2) and after an acute bout
at the end of it (T3). They were labelled High or Low Responder by a median split of a composite
hypertrophy score, eight each. 48 MS runs; three are dropped as outliers, leaving 45.

The primary comparison is the training interaction: whether the proteome changed differently over
training in HR than in LR. The study also asks which proteins, pathways and co-expression modules
separate the groups, and which track the phenotype measured on each biopsy.

## Every level runs contrasts, classify and associate

| Stage | Runs | Writes |
|---|---|---|
| [`00_Input/`](00_Input/README.md) | study data, no code | |
| [`01_Preprocess/`](01_Preprocess/README.md) | filtering, cyclic loess, limpa quantification | `proteins.rds`, `DAList_limpa.rds`, `DAList_corrected.rds` |
| [`02_Differential_Expression/`](02_Differential_Expression/README.md) | the protein level | `design.rds`, `fit.rds`, `associate.rds` |
| [`03_Pathway_Enrichment/`](03_Pathway_Enrichment/README.md) | the pathway level: gene sets, singscore | `gene_sets.rds`, `singscore.rds`, `set_tests.rds` |
| [`04_Network/`](04_Network/README.md) | the module level: WGCNA eigengenes | `modules.rds` |
| [`05_Summary/`](05_Summary/README.md) | permutation chance line, candidate tables across the three levels, exploratory panels | `01_permutation.xlsx`, `02_candidates.xlsx`, `03_panel.xlsx` |

Stages 02 to 04 share one layout. Step `02_Contrasts` tests the nine contrasts and a
baseline-adjusted ANCOVA, `03_Classify` scores eight classification tasks by AUC, and
`04_Associate` ties each feature to phenotype by a sample-level limma model, by change-score
correlation pooled and within each arm, by a test of whether that slope differs between arms, and
by a model with growth and strength side by side. Each level adds its own steps around them.
Each step's PDF opens with overview heatmaps of every feature, then gives every nominal result a
panel, and every table reports the tested, p < 0.05 and BH < 0.05 counts.

Each sub-stage holds `a_script/` (one R script), `b_reports/` (`<step>_figures.pdf`, 11 × 8.5 in),
`c_data/` (one workbook, plus the `.rds` a later step reads) and a README. Every workbook opens on a
`read_me` sheet and ends with `input_manifest` and `package_versions`. Data passes through disk, so
any sub-stage re-runs on its own once its inputs exist.

Every script sources `R/helpers.R`, which holds the statistics written for this study (blood
adjustment, ANCOVA, the Wilcoxon and AUC row, the Spearman p-values, the per-subject limma fit,
leave one subject out, Storey's q, the panel) and the heatmap, page and workbook plumbing the steps
share.

## One command per step, in order

```sh
Rscript -e 'renv::restore()'

for s in 01_Filtering/a_script/01_filter 02_Normalization/a_script/02_normalize \
         03_Quantification/a_script/03_quantify; do
  Rscript 01_Preprocess/$s.R
done

for s in 01_Design/a_script/01_design 02_Contrasts/a_script/02_contrasts \
         03_Classify/a_script/03_classify 04_Associate/a_script/04_associate; do
  Rscript 02_Differential_Expression/$s.R
done

for s in 00_Gene_Sets/a_script/00_gene_sets 01_Scores/a_script/01_scores \
         02_Contrasts/a_script/02_contrasts 03_Classify/a_script/03_classify \
         04_Associate/a_script/04_associate 05_Volcano/a_script/05_volcano \
         06_NES_Scatter/a_script/06_nes_scatter; do
  Rscript 03_Pathway_Enrichment/$s.R
done

for s in 01_Build_Modules/a_script/01_build_modules 02_Contrasts/a_script/02_contrasts \
         03_Classify/a_script/03_classify 04_Associate/a_script/04_associate \
         05_Characterise/a_script/05_characterise 06_Preserve/a_script/06_preserve; do
  Rscript 04_Network/$s.R
done

for s in 01_Permutation/a_script/01_permutation 02_Candidates/a_script/02_candidates \
         03_Panel/a_script/03_panel; do
  Rscript 05_Summary/$s.R
done
```

The run takes about 75 minutes on eight cores, most of it the panel, the preservation calibration
and the candidate permutation. `04_Network/05_Characterise` needs the STRING files; `00_Input/README.md` has the
download command.

## limpa models the missing values; blood is a covariate; BH never pools across screens

Filtering and normalisation use proteoDA (`DAList`, `zero_to_missing`, `filter_proteins_by_group`,
`filter_samples`, `normalize_data("cycloess")`). A missing value here says the protein was probably
faint, so limpa (Li, Cobbold & Smyth 2025) fits a detection probability curve and gives every
protein a value in every sample with a standard error (`dpcQuantByRow`, the vignette's route for
protein-level input). The protein fit is limpa's `dpcDE`: precision weights from those errors,
six cell means plus the blood index, subject blocked, robust moderation. Without precursor data
the gain over other methods is smaller than in a precursor pipeline such as BFR_limpa.

Blood enters every limma model as a covariate (index: mean log2 HBB and HBA1 against the sample
median). The rank tests, correlations, ANCOVA, singscore and WGCNA cannot carry a covariate, so they
read a copy with each protein's blood slope removed (`limma::removeBatchEffect`); the tests among
them use observed values only.

BH runs within each contrast, task, trait and term, or window and outcome. Storey's q sits beside
every BH column as a second exploratory correction over the same screen (`storey_q()` in
`R/helpers.R`: bootstrap pi0, equal to BH in screens under 100 tests). Both rank and flag nominal
results; neither gates them. Set tests pool the five
collections within a contrast; set classification and association split by collection.

## renv.lock pins every package

Stages 01 and 02: `proteoDA`, `limma`, `limpa`, `lme4`, `here`, `dplyr`,
`tidyr`, `tibble`, `purrr`, `stringr`, `forcats`, `readr`, `readxl`, `ggplot2`, `patchwork`, `pROC`,
`qvalue`, `ComplexHeatmap`, `circlize`, `writexl`, `sessioninfo`. Stage 03 adds `fgsea`, `singscore`, `msigdbr`, `GO.db`,
`GSEABase`, `AnnotationDbi`, `ggrepel` and `enrichVolcano`. Stage 04 adds `WGCNA`,
`clusterProfiler`, `enrichplot`, `STRINGdb`, `tidygraph`, `ggraph` and `igraph`. Stage 05 adds
`nestedcv`, which brings `glmnet` and `caret`.

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

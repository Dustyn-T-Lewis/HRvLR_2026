# HRvLR_F1_setup.R — Shared setup for Figure 1 ----
# Provides: norm_df, imp_df, dep_df, meta, samp_names (norm), imp_samps (imp),
#           ann_cols, imp_mat, BEST_IMP_METHOD, MAIN_CONTRASTS, RPT_DIR, DAT_DIR
# Plus all style.R exports (palettes, themes, helpers)

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ComplexHeatmap)
  library(fgsea)
  library(msigdbr)
  library(grid)
  library(vegan)
})

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

# Paths
NORM_FILE <- "01_normalization/c_data/02_normalized.csv"
IMP_FILE  <- "02_Imputation/c_data/01_imputed.csv"
DEP_FILE  <- "03_DEP/c_data/03_combined_results.csv"
META_FILE <- "00_input/HRvLR_meta.csv"

RPT_DIR <- "04_Figures/F1/b_reports"
DAT_DIR <- "04_Figures/F1/c_data"
dir.create(RPT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(RPT_DIR, "supplementary"), showWarnings = FALSE)

MAIN_CONTRASTS <- c("Training_HR", "Training_LR", "Acute_HR",
                    "Acute_LR", "Baseline_HRvLR")

DISPLAY_ORDER <- c("Baseline_HRvLR",
                   "Training_HR", "Training_LR", "Training_Interaction",
                   "Acute_HR", "Acute_LR", "Acute_Interaction")

# Load data
norm_df <- read_csv(NORM_FILE, show_col_types = FALSE)
imp_df  <- read_csv(IMP_FILE,  show_col_types = FALSE)
dep_df  <- read_csv(DEP_FILE,  show_col_types = FALSE)

ann_cols   <- c("uniprot_id", "protein", "gene", "description")
samp_names <- setdiff(names(norm_df), ann_cols)

# Imputed data may have fewer samples (QC-removed during imputation)
imp_ann    <- intersect(names(imp_df), ann_cols)
imp_samps  <- setdiff(names(imp_df), imp_ann)

# Metadata — from CSV, joined to imp_samps (the analysis-ready sample set)
meta_raw <- read_csv(META_FILE, show_col_types = FALSE)
meta <- tibble(sample_id = imp_samps) |>
  left_join(meta_raw |> select(Col_ID, Subject_ID, Group, Timepoint, Group_Time),
            by = c("sample_id" = "Col_ID")) |>
  mutate(
    subject   = Subject_ID,
    group     = factor(Group_Time,
                       levels = c("HR_T1", "HR_T2", "HR_T3",
                                  "LR_T1", "LR_T2", "LR_T3")),
    Group     = factor(Group, levels = c("HR", "LR")),
    Timepoint = factor(Timepoint, levels = c("T1", "T2", "T3"))
  )

cat(sprintf("Loaded: %d norm proteins (%d samples), %d imp proteins (%d samples), %d DEP rows\n",
            nrow(norm_df), length(samp_names), nrow(imp_df), length(imp_samps), nrow(dep_df)))

# Best imputation method
imp_summary <- readLines("02_Imputation/c_data/09_imputation_summary.txt")
BEST_IMP_METHOD <- toupper(trimws(sub(".*=\\s*", "",
  grep("^best_method", imp_summary, value = TRUE))))
if (length(BEST_IMP_METHOD) == 0) BEST_IMP_METHOD <- "IMPUTED"

# Imputed matrix (proteins × samples)
imp_mat <- as.matrix(imp_df[, imp_samps])
rownames(imp_mat) <- imp_df$gene

# Cliff's delta (non-parametric effect size) — used by panels A and B
cliffs_delta <- function(x, y) {
  d <- outer(x, y, function(a, b) sign(a - b))
  sum(d) / (length(x) * length(y))
}

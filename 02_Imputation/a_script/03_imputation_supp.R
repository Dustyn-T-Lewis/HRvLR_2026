#!/usr/bin/env Rscript
# Rebuild the HRvLR Stage 02 supplementary workbook from top-level outputs.

library(openxlsx)
library(readr)
library(tibble)

setwd(rprojroot::find_rstudio_root_file())

DATA_DIR <- "02_Imputation/c_data"
BENCH_DIR <- file.path(DATA_DIR, "benchmark")
OUT_FILE <- file.path(DATA_DIR, "02_imputation.xlsx")

miss_class <- read_csv(file.path(DATA_DIR, "02_mar_mnar_classification.csv"),
  show_col_types = FALSE)
imputed_df <- read_csv(file.path(DATA_DIR, "01_imputed.csv"),
  show_col_types = FALSE)
mask_df <- read_csv(file.path(DATA_DIR, "07_imputation_mask.csv"),
  show_col_types = FALSE)
mnar_audit <- read_csv(file.path(DATA_DIR, "08_mnar_imputation_audit.csv"),
  show_col_types = FALSE)
summary_lines <- readLines(file.path(DATA_DIR, "09_imputation_summary.txt"))

summary_df <- tibble(
  metric = sub(" = .*", "", summary_lines),
  value = sub(".* = ", "", summary_lines)
)

readme_df <- tibble(
  sheet = c(
    "imputed_matrix",
    "mar_mnar_classification",
    "imputation_mask",
    "mnar_audit",
    "summary",
    "benchmark_ranking",
    "benchmark_summary"
  ),
  description = c(
    "Fully imputed matrix written by the canonical YvO-style missForest workflow.",
    "Per-feature MAR/MNAR calls, vote breakdown, and reliability flags.",
    "Boolean mask marking values that were originally missing.",
    "Pre/post mean shift audit for MNAR-classified features.",
    "Key pipeline parameters and output counts.",
    "Optional benchmark composite ranking if exploratory benchmark outputs exist.",
    "Optional benchmark NRMSE/PSS summary if exploratory benchmark outputs exist."
  )
)

add_sheet <- function(wb, name, df) {
  header_style <- createStyle(
    textDecoration = "bold",
    border = "Bottom",
    fgFill = "#DCE6F1"
  )
  addWorksheet(wb, name)
  writeData(wb, name, df, headerStyle = header_style)
  freezePane(wb, name, firstActiveRow = 2)
  setColWidths(wb, name, cols = seq_len(ncol(df)), widths = "auto")
}

wb <- createWorkbook()
add_sheet(wb, "README", readme_df)
add_sheet(wb, "imputed_matrix", imputed_df)
add_sheet(wb, "mar_mnar_classification", miss_class)
add_sheet(wb, "imputation_mask", mask_df)
add_sheet(wb, "mnar_audit", mnar_audit)
add_sheet(wb, "summary", summary_df)

bench_rank_path <- file.path(BENCH_DIR, "04_composite_ranking.csv")
if (file.exists(bench_rank_path)) {
  add_sheet(
    wb,
    "benchmark_ranking",
    read_csv(bench_rank_path, show_col_types = FALSE)
  )
}

bench_sum_path <- file.path(BENCH_DIR, "03_benchmark_summary.csv")
if (file.exists(bench_sum_path)) {
  add_sheet(
    wb,
    "benchmark_summary",
    read_csv(bench_sum_path, show_col_types = FALSE)
  )
}

saveWorkbook(wb, OUT_FILE, overwrite = TRUE)
cat(sprintf("Supplementary workbook written: %s\n", OUT_FILE))

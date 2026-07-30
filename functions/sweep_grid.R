# Screen axes and shared leaf plumbing for the exploratory sweep. One place
# defines the feature levels, timepoint configs, and method sets so every root
# orchestrator and composite reads the same grid.

pacman::p_load(here, dplyr, openxlsx)

SWEEP_LEVELS <- c("pathways", "modules", "proteins")
SWEEP_LEVEL_KEY <- c(
  pathways = "singscore", modules = "eigengenes", proteins = "proteins"
)
SWEEP_LEVEL_LABEL <- c(
  pathways = "pathways (singscore)", modules = "modules (ME)*",
  proteins = "proteins*"
)

SWEEP_CONFIGS <- c(
  "T1", "T2", "T3", "training", "acute", "total", "trajectory"
)
CONFIG_ROLE <- c(
  T1 = "baseline forecast", T2 = "separability", T3 = "separability",
  training = "training change", acute = "acute change",
  total = "whole-study change", trajectory = "concatenated"
)

ADAPT_OUTCOMES <- c(
  "comp_hypertrophy", "d_fcsa_I", "d_fcsa_II", "d_mcsa",
  "d_1rm_legpress", "d_1rm_ext"
)

ASSOC_CONT_METHODS <- c("limma", "lm", "spearman")
ASSOC_GROUP_METHODS <- c("limma", "wilcoxon")
CLASS_METHODS <- c("enet", "lasso", "ridge", "spls", "pam", "rf", "svm")
CONT_METHODS <- c("enet", "lasso", "ridge", "spls", "rf", "svm")
PLAIN_METHOD <- "plain"

# Methods that run at a given level: the unpenalised plain model is only valid
# on the low-dimension module space (p < n).
methods_for_level <- function(base_methods, level) {
  if (level == "modules") c(base_methods, PLAIN_METHOD) else base_methods
}

sweep_root_dir <- function(root) here("05_Figures", root)

sweep_leaf_dir <- function(root, level, config, method) {
  d <- file.path(sweep_root_dir(root), level, config, method)
  for (sub in c("b_reports", "c_data")) {
    dir.create(file.path(d, sub), recursive = TRUE, showWarnings = FALSE)
  }
  d
}

# A leaf is done when its workbook exists, so a killed run resumes by skipping
# the leaves on disk.
leaf_done <- function(root, level, config, method) {
  file.exists(file.path(
    sweep_root_dir(root), level, config, method, "c_data", "results.xlsx"
  ))
}

write_sweep_workbook <- function(path, sheets) {
  wb <- createWorkbook()
  for (nm in names(sheets)) {
    addWorksheet(wb, nm)
    writeData(wb, nm, sheets[[nm]])
  }
  saveWorkbook(wb, path, overwrite = TRUE)
}

# A lead clears its permutation null AND beats the trivial baseline: predicting
# the mean for Q2, chance for AUC. Ridge and the unpenalised plain model
# collapse their nulls, so p alone promotes cells whose metric is worse than
# the baseline. Every place that counts or highlights a lead calls this, so the
# manifest, the roll-up, the composites and the specification curve cannot
# disagree. Association is in-sample with no permutation null; it returns NA.
LEAD_BASELINE <- c(cont = 0, class = 0.5)

# Models that retain every coefficient, so their fold-selection frequencies
# describe the whole feature space rather than a signature. Both leaf-panel
# builders read this; with two copies they disagreed, and the pre-split panel
# was ranking ridge features alphabetically under the label "top features".
DENSE_MODELS <- "ridge"

# A cell is one level x config x model, and one outcome where the root sweeps
# several. Classification workbooks carry no outcome column, so the key is
# whichever of these the table actually has.
cell_key <- function(df) {
  intersect(c("level", "config", "outcome", "model"), names(df))
}

# Each cell reports at its own best resolution. Taking max(B) over the pooled
# table instead would keep only the cells that reached the highest B anywhere,
# silently dropping every cell swept at a lower B.
best_b_per_cell <- function(df) {
  slice_max(df, .data$B, by = all_of(cell_key(df)), with_ties = FALSE)
}

is_lead_at <- function(metric, p, baseline) {
  !is.na(p) & p < 0.05 & metric > baseline
}

is_lead <- function(metric, p, kind) {
  if (!kind %in% names(LEAD_BASELINE)) {
    return(rep(NA, length(p)))
  }
  is_lead_at(metric, p, LEAD_BASELINE[[kind]])
}

# Which metric and p a root's summary carries, so callers need not hardcode it.
root_kind <- function(root) {
  if (grepl("classification", root)) "class" else "cont"
}

root_metric_col <- function(kind) if (kind == "class") "estimate" else "q2"

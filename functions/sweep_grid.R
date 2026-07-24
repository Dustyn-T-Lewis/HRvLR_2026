# Screen axes and shared leaf plumbing for the exploratory sweep. One place
# defines the feature levels, timepoint configs, and method sets so every root
# orchestrator and composite reads the same grid.

pacman::p_load(here, openxlsx)

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

sweep_root_dir <- function(root) here("04_Figures", root)

sweep_leaf_dir <- function(root, level, config, method) {
  d <- file.path(sweep_root_dir(root), level, config, method)
  for (sub in c("a_script", "b_reports", "c_data")) {
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

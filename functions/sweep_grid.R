# Screen axes and shared leaf plumbing for the exploratory sweep. One place
# defines the feature levels, timepoint configs, and method sets so every root
# orchestrator and composite reads the same grid.

pacman::p_load(here, dplyr, readr, openxlsx, digest)

SWEEP_LEVELS <- c("pathways", "modules", "proteins")
SWEEP_LEVEL_KEY <- c(
  pathways = "singscore", modules = "eigengenes", proteins = "proteins"
)
SWEEP_LEVEL_LABEL <- c(
  pathways = "pathways (singscore)", modules = "modules (ME)*",
  proteins = "proteins*"
)
# The level names are the vocabulary of the screen; the directories they answer
# in are not. Figures that read a level's results table need the mapping, and
# renaming the levels themselves would change factor levels that reach the plots.
SWEEP_LEVEL_DIR <- c(
  pathways = "02_Pathways", modules = "03_WGCNA", proteins = "01_Proteins"
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

CLASS_METHODS <- c("enet", "lasso", "ridge", "spls", "pam", "rf", "svm")
CONT_METHODS <- c("enet", "lasso", "ridge", "spls", "rf", "svm")
PLAIN_METHOD <- "plain"

# Methods that run at a given level: the unpenalised plain model is only valid
# on the low-dimension module space (p < n).
methods_for_level <- function(base_methods, level) {
  if (level == "modules") c(base_methods, PLAIN_METHOD) else base_methods
}

sweep_root_dir <- function(root) here("04_Figures", root)

# Every cell in a root lives in four tables under <root>/c_data, one per sheet,
# each keyed by (level, config, phenotype, model). This replaces the directory
# per cell and the second copy the split step used to make.
#
# One table per sheet rather than one table for the root, because the four
# sheets are four different grains: a summary row per B, a null row per
# permutation draw, a prediction row per subject, a selection row per feature.
# Unioning them pads every row with the other three sheets' columns and puts
# the observed Q2 and its permutation draws in one `q2` column separated only
# by a discriminator. Four dense tables need no padding and no discriminator.
SWEEP_SHEETS <- c("summary", "null", "predictions", "selection")

# Added on write and stripped on read, so a cell reads back as the sheet it was
# written from. `model` is not here: all four sheets carry it already.
SWEEP_STORE_COLS <- c("level", "config", "phenotype", "fingerprint")

sweep_store_path <- function(root, sheet, root_dir = sweep_root_dir(root)) {
  file.path(root_dir, "c_data", paste0("cells_", sheet, ".csv"))
}

read_sweep_store <- function(root, sheet, root_dir = sweep_root_dir(root)) {
  path <- sweep_store_path(root, sheet, root_dir)
  if (!file.exists(path)) {
    return(NULL)
  }
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}

# A sheet is filtered only when it carries the outcome key; the classification
# null, predictions and selection sheets do not, so they copy across whole.
filter_to_outcome <- function(df, outcome) {
  if (is.null(df) || !"outcome" %in% names(df)) {
    return(df)
  }
  df[df$outcome == outcome, , drop = FALSE]
}

# Cell order is load-bearing: the composites rank cells and take the top 12, so
# ties are broken by whatever order the table arrives in. Sorting on write makes
# that independent of the order the sweep computed cells in.
#
# `method = "radix"` sorts in the C locale, which is the whole point. The order
# used to come from Sys.glob() over the leaf directories, so it inherited the
# session's collation -- en_US puts `acute` ahead of `T1`, the C locale does the
# reverse -- and the same code on a Linux box could therefore pick a different
# top 12. Radix is also stable, so rows within a cell keep the order they were
# computed in.
sort_sweep_store <- function(rows) {
  ord <- order(
    rows$level, rows$config, rows$phenotype, rows$model,
    method = "radix"
  )
  rows[ord, , drop = FALSE]
}

# TRUE for the store rows belonging to one cell.
cell_rows_of <- function(store, level, config, phenotype, model) {
  store$level == level & store$config == config &
    store$phenotype == phenotype & store$model == model
}

# Rewriting the cell's rows rather than appending keeps a refit from stacking a
# second copy on top of the first.
write_sweep_cell <- function(root, level, config, phenotype, model, sheets,
                             fingerprint, root_dir = sweep_root_dir(root)) {
  dir.create(file.path(root_dir, "c_data"),
    recursive = TRUE, showWarnings = FALSE
  )
  for (nm in names(sheets)) {
    rows <- as.data.frame(sheets[[nm]]) |>
      dplyr::mutate(
        level = level, config = config, phenotype = phenotype,
        model = model, fingerprint = fingerprint
      )
    store <- read_sweep_store(root, nm, root_dir)
    if (!is.null(store)) {
      keep <- !cell_rows_of(store, level, config, phenotype, model)
      rows <- dplyr::bind_rows(store[keep, , drop = FALSE], rows)
    }
    readr::write_csv(
      sort_sweep_store(rows), sweep_store_path(root, nm, root_dir)
    )
  }
  invisible(TRUE)
}

read_sweep_cell <- function(root, level, config, phenotype, model, sheet,
                            root_dir = sweep_root_dir(root)) {
  store <- read_sweep_store(root, sheet, root_dir)
  if (is.null(store)) {
    return(NULL)
  }
  hit <- store[cell_rows_of(store, level, config, phenotype, model), ,
    drop = FALSE
  ]
  dplyr::select(hit, -dplyr::any_of(SWEEP_STORE_COLS))
}

# A leaf is done when the store holds its rows AND they carry the fingerprint
# of the input it was fitted on, so a killed run still resumes from disk while
# a change upstream forces a refit. Presence alone was not enough: when the
# protein set moved on 2026-07-30 every one of the 945 leaves reported
# "skip (done)" against results fitted on a matrix that no longer existed, and
# nothing said so.
#
# The permutation grid is part of the input. A leaf swept at B = 0 carries only
# point estimates, and without b_grid here the later B = 200 pass would take it
# as done and the sweep would finish with no nulls at all.
# Any analysis that caches a result needs to record what it was computed from.
input_fingerprint <- function(...) {
  substr(digest::digest(list(...)), 1, 12)
}

sweep_fingerprint <- function(bundle, b_grid) {
  input_fingerprint(bundle$feature_sets, as.integer(b_grid))
}

leaf_done <- function(root, level, config, method, fingerprint,
                      root_dir = sweep_root_dir(root)) {
  store <- read_sweep_store(root, "summary", root_dir)
  if (is.null(store)) {
    return(FALSE)
  }
  hit <- dplyr::filter(
    store,
    .data$level == !!level, .data$config == !!config,
    .data$model == !!method
  )
  nrow(hit) > 0L && all(hit$fingerprint == fingerprint)
}

write_sweep_workbook <- function(path, sheets, fingerprint = NULL) {
  if (!is.null(fingerprint)) {
    sheets$provenance <- data.frame(fingerprint = fingerprint)
  }
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

# Each cell reports at its own best B. Filtering on max(B) over the pooled
# table instead would silently drop every cell swept at a lower B -- the
# documented split runs the fast levels at 0/200/1000 and proteins at 0/200,
# so that filter deletes the entire protein level from the figure with no gap
# and no warning.
root_cells <- function(root) {
  best_b_per_cell(read_sweep_store(root, "summary"))
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

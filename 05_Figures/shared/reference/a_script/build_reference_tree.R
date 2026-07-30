# A worked reference for every cell type in the screen: one panel per
# stage x feature level x timepoint config, written into a tree that mirrors
# F05/F06. Each reference is drawn from that cell's real output, so it shows
# both the intended layout and what the data actually supports there.
#
# Stage sets the statistic (classification: AUC; prediction: Q^2), the feature
# level sets how a driver is named (gene
# symbol, pathway name, module ORA term), and the config sets the readout
# label (abundance, logFC, trajectory).

suppressWarnings(suppressMessages({
  library(here)
  source(here("functions", "sweep_composites.R"))
  source(here("functions", "sweep_pred_leaf.R"))
  source(here("functions", "sweep_drivers.R"))
  source(here("functions", "sweep_cell_panel.R"))
  pacman::p_load(openxlsx, dplyr, ggplot2, patchwork, stringr, scales)
}))

REF_ROOT <- here("05_Figures", "shared", "reference")
SPARSE_MODELS <- c("lasso", "enet", "spls", "pam")
LEAD_OUTCOME <- "d_mcsa"

STAGES <- c(
  F05_classification = "classification",
  F06_prediction = "prediction"
)

# T1/T2/T3 are levels, the paired configs are changes, trajectory concatenates.
readout_of <- function(config, level) {
  base <- switch(config,
    T1 = ,
    T2 = ,
    T3 = "abundance",
    training = ,
    acute = ,
    total = "logFC",
    trajectory = "trajectory"
  )
  switch(level,
    pathways = paste(base, "(pathway activity)"),
    modules = paste(base, "(module eigengene / kME)"),
    base
  )
}

base_feature <- function(x) sub("@T[123]$", "", x)

module_ora_label <- function() {
  read.xlsx(
    here("04_Features", "modules", "c_data", "WGCNA_source_data.xlsx"),
    "module_ora"
  ) |>
    group_by(.data$module) |>
    slice_min(.data$padj, n = 1, with_ties = FALSE) |>
    ungroup() |>
    transmute(
      feature = paste0("ME_", .data$module),
      label = ifelse(
        is.na(.data$term), sprintf("%s: no enriched term", .data$module),
        sprintf("%s: %s", .data$module, str_trunc(.data$term, 30))
      )
    )
}

ORA_LABELS <- module_ora_label()

feature_label <- function(features, level) {
  if (level == "modules") {
    out <- ORA_LABELS$label[match(features, ORA_LABELS$feature)]
    ifelse(is.na(out), features, out)
  } else if (level == "pathways") {
    clean_pathway_name(features)
  } else {
    features
  }
}

cell_files <- function(stage, level, config) {
  Sys.glob(file.path(
    sweep_root_dir(stage), level, config, "*", "*", "c_data", "results.xlsx"
  ))
}

# A leaf's phenotype directory sits two levels above its method directory, so
# a file already matched by `cell_files()` can be filtered down to one outcome
# without re-globbing.
file_phenotype <- function(files) {
  vapply(
    files, function(f) basename(dirname(dirname(dirname(f)))), character(1)
  )
}

# Prediction / classification drivers: fold-selection frequency pooled over the
# sparse models only, since ridge and the dense learners carry no signature.
pred_drivers <- function(stage, level, config, n = 10) {
  rows <- lapply(cell_files(stage, level, config), function(f) {
    sel <- read.xlsx(f, "selection")
    if (!nrow(sel)) {
      return(NULL)
    }
    if ("outcome" %in% names(sel)) {
      sel <- filter(sel, .data$outcome == LEAD_OUTCOME)
    }
    filter(sel, .data$model %in% SPARSE_MODELS)
  })
  out <- bind_rows(rows)
  if (!nrow(out)) {
    return(out)
  }
  out |>
    mutate(feature = base_feature(.data$feature)) |>
    group_by(.data$feature) |>
    summarise(
      cells = dplyr::n(), score = mean(.data$freq), .groups = "drop"
    ) |>
    slice_max(.data$score, n = n, with_ties = FALSE) |>
    mutate(label = feature_label(.data$feature, level))
}

empty_panel <- function(msg) {
  ggplot() +
    annotate("text", 0, 0,
      label = msg, size = 3, colour = "grey45",
      fontface = "italic"
    ) +
    theme_void()
}

driver_panel <- function(d, level, xlab, subtitle, signed = FALSE) {
  driver_bars(d, level, xlab, "Named drivers", subtitle, signed = signed)
}

# Classification reads best as arm separation: the out-of-fold score each
# subject received, split by the arm they actually belong to. Overlap is the
# null; a gap is the signal.
arm_separation_panel <- function(stage, level, config) {
  files <- cell_files(stage, level, config)
  summ <- bind_rows(lapply(files, function(f) {
    cbind(file = f, read.xlsx(f, "summary"))
  })) |>
    best_b_per_cell()
  top <- summ |> slice_min(.data$perm_p, n = 1, with_ties = FALSE)
  pr <- read.xlsx(top$file, "predictions") |>
    mutate(arm = ifelse(.data$y == 1, "HR", "LR"))

  ggplot(pr, aes(.data$arm, .data$pred, fill = .data$arm)) +
    geom_boxplot(
      width = 0.5, alpha = 0.35, outlier.shape = NA,
      colour = "grey30"
    ) +
    geom_jitter(aes(colour = .data$arm), width = 0.12, size = 2, alpha = 0.9) +
    scale_fill_manual(values = GROUP_COLORS, guide = "none") +
    scale_colour_manual(values = GROUP_COLORS, guide = "none") +
    annotate("label",
      x = 1.5, y = Inf, vjust = 1.1,
      label = cell_footer("AUC", top$estimate, top$perm_p, stage, level),
      size = 2.6, fill = alpha("white", 0.9), lineheight = 0.95
    ) +
    labs(
      x = NULL, y = "out-of-fold score",
      title = sprintf("Arm separation -- %s", top$model),
      subtitle = "Held-out score by true arm. Overlap = null."
    ) +
    FIG_THEME +
    theme(plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE))
}

# The statistics block: the stage's metric against its null, with the raw
# permutation p taken within the cell.
stat_panel <- function(stage, level, config) {
  files <- cell_files(stage, level, config)
  if (stage == "F05_classification") {
    return(arm_separation_panel(stage, level, config))
  }

  is_class <- FALSE
  summ <- bind_rows(lapply(files, function(f) read.xlsx(f, "summary"))) |>
    best_b_per_cell()
  if (!is_class) summ <- filter(summ, .data$outcome == LEAD_OUTCOME)
  metric <- if (is_class) "estimate" else "q2"
  pcol <- if (is_class) "perm_p" else "perm_p_q2"
  top <- summ |> slice_min(.data[[pcol]], n = 1, with_ties = FALSE)

  null <- bind_rows(lapply(files, function(f) {
    nl <- read.xlsx(f, "null")
    if ("outcome" %in% names(nl)) filter(nl, .data$outcome == LEAD_OUTCOME) else nl
  }))
  null <- filter(null, .data$model == top$model)
  vcol <- if (is_class) "auc" else "q2"

  ggplot(null, aes(.data[[vcol]])) +
    geom_histogram(
      bins = 28, fill = "grey78", colour = "white",
      linewidth = 0.1
    ) +
    geom_vline(
      xintercept = if (is_class) 0.5 else 0, linetype = "dotted",
      colour = "grey45"
    ) +
    geom_vline(xintercept = top[[metric]], colour = "#B2182B", linewidth = 1) +
    annotate("label",
      x = -Inf, y = Inf, hjust = -0.05, vjust = 1.15,
      label = cell_footer(
        if (is_class) "AUC" else "Q2", top[[metric]], top[[pcol]], stage, level
      ),
      size = 2.6, fill = alpha("white", 0.9), lineheight = 0.95
    ) +
    labs(
      x = sprintf("permutation null (%s)", if (is_class) "AUC" else "Q2"),
      y = "count", title = sprintf("Out-of-sample statistic -- %s", top$model),
      subtitle = "Observed (red) vs its own null."
    ) +
    FIG_THEME +
    theme(plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE))
}

build_reference <- function(stage, level, config) {
  if (!length(cell_files(stage, level, config))) {
    return(invisible(NULL))
  }
  d <- pred_drivers(stage, level, config)
  target <- if (stage == "F05_classification") "HR vs LR" else LEAD_OUTCOME
  sub <- sprintf("Selected for %s across sparse models", target)

  panel <- (driver_panel(d, level, "mean fold-selection freq.", sub) |
    stat_panel(stage, level, config)) +
    plot_annotation(
      title = sprintf(
        "REFERENCE -- %s -- %s -- %s", STAGES[[stage]], level, config
      ),
      subtitle = sprintf(
        "Readout: %s. %s Drivers carry real names; the statistic carries its null.",
        readout_of(config, level), LEAK_NOTE[[level]]
      ),
      tag_levels = "A",
      theme = theme(
        plot.title = element_text(face = "bold", size = FIG_TITLE_SIZE + 1),
        plot.subtitle = element_text(
          face = "italic", size = FIG_SUBTITLE_SIZE, colour = "grey30"
        ),
        plot.tag = element_text(face = "bold", size = FIG_TAG_SIZE)
      )
    )

  out <- file.path(REF_ROOT, stage, level, config)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  save_panel(panel, file.path(out, "reference"), width = 230, height = 95)
  invisible(TRUE)
}

grid <- expand.grid(
  stage = names(STAGES), level = SWEEP_LEVELS, config = SWEEP_CONFIGS,
  stringsAsFactors = FALSE
)
n_built <- sum(vapply(seq_len(nrow(grid)), function(i) {
  g <- grid[i, ]
  ok <- !is.null(build_reference(g$stage, g$level, g$config))
  if (ok) message(sprintf("[ref] %s | %s -- %s", g$stage, g$level, g$config))
  ok
}, logical(1)))

message(sprintf("built %d cell references under %s", n_built, REF_ROOT))

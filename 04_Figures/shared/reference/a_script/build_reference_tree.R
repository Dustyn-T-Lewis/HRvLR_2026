# A worked reference for every cell type in the screen: one panel per
# stage x feature level x timepoint config, written into a tree that mirrors
# F04/F05/F06. Each reference is drawn from that cell's real output, so it shows
# both the intended layout and what the data actually supports there.
#
# Stage sets the statistic (association: effect and BH q; classification: AUC;
# prediction: Q^2), the feature level sets how a driver is named (gene symbol,
# pathway name, module ORA term), and the config sets the readout label
# (abundance, logFC, trajectory).

suppressWarnings(suppressMessages({
  library(here)
  source(here("functions", "sweep_composites.R"))
  source(here("functions", "sweep_pred_leaf.R"))
  pacman::p_load(openxlsx, dplyr, ggplot2, patchwork, stringr, scales)
}))

REF_ROOT <- here("04_Figures", "shared", "reference")
SPARSE_MODELS <- c("lasso", "enet", "spls", "pam")
LEAD_OUTCOME <- "d_mcsa"

STAGES <- c(
  F04_association = "association",
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
    here("04_Figures", "shared", "WGCNA", "c_data", "WGCNA_source_data.xlsx"),
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
    sweep_root_dir(stage), level, config, "*", "c_data", "results.xlsx"
  ))
}

# Association drivers: strongest features by raw p in the lead outcome, pooled
# over the methods run in that cell.
assoc_drivers <- function(stage, level, config, n = 10) {
  rows <- lapply(cell_files(stage, level, config), function(f) {
    if (!LEAD_OUTCOME %in% getSheetNames(f)) {
      return(NULL)
    }
    read.xlsx(f, LEAD_OUTCOME) |>
      mutate(method = basename(dirname(dirname(f))))
  })
  bind_rows(rows) |>
    mutate(feature = base_feature(.data$feature)) |>
    group_by(.data$feature) |>
    summarise(
      score = -log10(min(.data$p, na.rm = TRUE)),
      best_bh = min(.data$bh, na.rm = TRUE), .groups = "drop"
    ) |>
    slice_max(.data$score, n = n, with_ties = FALSE) |>
    mutate(label = feature_label(.data$feature, level))
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

driver_panel <- function(d, level, xlab, subtitle) {
  if (!nrow(d)) {
    return(empty_panel("no sparse signature\nin this cell"))
  }
  ggplot(d, aes(.data$score, stats::reorder(.data$label, .data$score))) +
    geom_col(fill = SPEC_LEVEL_COLORS[[level]], width = 0.72, alpha = 0.9) +
    labs(x = xlab, y = NULL, title = "Named drivers", subtitle = subtitle) +
    FIG_THEME +
    theme(
      axis.text.y = element_text(size = 6.5),
      plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE)
    )
}

# The statistics block: the stage's metric against its null, with p and the BH q
# taken within the cell; bold when q clears 0.05.
stat_panel <- function(stage, level, config) {
  files <- cell_files(stage, level, config)
  summ <- bind_rows(lapply(files, function(f) read.xlsx(f, "summary")))
  if (stage == "F04_association") {
    d <- bind_rows(lapply(files, function(f) {
      if (!LEAD_OUTCOME %in% getSheetNames(f)) NULL else read.xlsx(f, LEAD_OUTCOME)
    }))
    n_bh <- sum(d$bh < 0.05, na.rm = TRUE)
    return(
      ggplot(d, aes(.data$effect, -log10(.data$p))) +
        geom_point(size = 0.7, alpha = 0.5, colour = SPEC_LEVEL_COLORS[[level]]) +
        labs(
          x = sprintf("effect vs %s", LEAD_OUTCOME),
          y = expression(-log[10] * " p"), title = "In-sample statistic",
          subtitle = sprintf("BH q<.05 in cell: %d / %d", n_bh, nrow(d))
        ) +
        FIG_THEME +
        theme(plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE))
    )
  }

  is_class <- stage == "F05_classification"
  summ <- filter(summ, .data$B == max(.data$B))
  if (!is_class) summ <- filter(summ, .data$outcome == LEAD_OUTCOME)
  metric <- if (is_class) "estimate" else "q2"
  pcol <- if (is_class) "perm_p" else "perm_p_q2"
  summ <- summ |> mutate(q = stats::p.adjust(.data[[pcol]], "BH"))
  top <- summ |> slice_min(.data[[pcol]], n = 1, with_ties = FALSE)

  null <- bind_rows(lapply(files, function(f) {
    nl <- read.xlsx(f, "null")
    if ("outcome" %in% names(nl)) filter(nl, .data$outcome == LEAD_OUTCOME) else nl
  }))
  null <- filter(null, .data$model == top$model)
  vcol <- if (is_class) "auc" else "q2"
  sig <- top$q < 0.05

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
      label = sprintf(
        "%s = %.2f\np = %.3f\nq = %.3f",
        if (is_class) "AUC" else "Q2", top[[metric]], top[[pcol]], top$q
      ),
      size = 2.6, fontface = if (sig) "bold" else "plain",
      fill = alpha("white", 0.9), lineheight = 0.95
    ) +
    labs(
      x = sprintf("permutation null (%s)", if (is_class) "AUC" else "Q2"),
      y = "count", title = sprintf("Out-of-sample statistic -- %s", top$model),
      subtitle = "Observed (red) vs its own null. Bold = q < .05."
    ) +
    FIG_THEME +
    theme(plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE))
}

build_reference <- function(stage, level, config) {
  if (!length(cell_files(stage, level, config))) {
    return(invisible(NULL))
  }
  is_assoc <- stage == "F04_association"
  d <- if (is_assoc) {
    assoc_drivers(stage, level, config)
  } else {
    pred_drivers(stage, level, config)
  }
  xlab <- if (is_assoc) "-log10 p (best method)" else "mean fold-selection freq."
  sub <- if (is_assoc) {
    sprintf("Strongest %s associations, named", LEAD_OUTCOME)
  } else {
    sprintf("Selected for %s across sparse models", LEAD_OUTCOME)
  }

  panel <- (driver_panel(d, level, xlab, sub) | stat_panel(stage, level, config)) +
    plot_annotation(
      title = sprintf(
        "REFERENCE -- %s | %s | %s", STAGES[[stage]], level, config
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
  if (ok) message(sprintf("[ref] %s | %s | %s", g$stage, g$level, g$config))
  ok
}, logical(1)))

message(sprintf("built %d cell references under %s", n_built, REF_ROOT))

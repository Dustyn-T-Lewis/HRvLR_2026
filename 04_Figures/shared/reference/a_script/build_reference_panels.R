# Design references for the figure review. Every panel is drawn from the real
# sweep output, so the layouts double as a check that the underlying data
# supports them. Three references: the main-composite idiom (named features and
# the per-cell statistics block), the per-figure validation supplement, and the
# shared upstream QC figure.

suppressWarnings(suppressMessages({
  library(here)
  source(here("functions", "sweep_composites.R"))
  source(here("functions", "sweep_pred_leaf.R"))
  source(here("functions", "pred_features.R"))
  source(here("functions", "sweep_cell_panel.R"))
  pacman::p_load(openxlsx, dplyr, tidyr, ggplot2, patchwork, stringr, scales)
}))

REF_DIR <- here("04_Figures", "shared", "reference", "b_reports")
dir.create(REF_DIR, recursive = TRUE, showWarnings = FALSE)

wgcna_book <- here(
  "04_Figures", "shared", "WGCNA", "c_data",
  "WGCNA_source_data.xlsx"
)

# Modules carry no biology in their label, so each is renamed by its strongest
# over-represented term. Proteins and pathways already self-name.
module_ora_label <- function() {
  read.xlsx(wgcna_book, "module_ora") |>
    group_by(.data$module) |>
    slice_min(.data$padj, n = 1, with_ties = FALSE) |>
    ungroup() |>
    transmute(
      feature = paste0("ME_", .data$module),
      label = ifelse(
        is.na(.data$term),
        sprintf("%s: no enriched term", .data$module),
        sprintf("%s: %s", .data$module, str_trunc(.data$term, 32))
      )
    )
}

# The trajectory encoding suffixes each feature with its timepoint, so pooling
# across configs must collapse back to the base feature: one pathway, not three.
base_feature <- function(x) sub("@T[123]$", "", x)

feature_label <- function(features, level) {
  if (level == "modules") {
    lab <- module_ora_label()
    out <- lab$label[match(features, lab$feature)]
    ifelse(is.na(out), features, out)
  } else if (level == "pathways") {
    clean_pathway_name(features)
  } else {
    features
  }
}

# Models that yield a genuine signature. ridge keeps every coefficient non-zero,
# so it "selects" the whole feature space and would swamp the pooled frequency
# with a flat 100%; it is excluded, as are the dense learners that record nothing.
SPARSE_MODELS <- c("lasso", "enet", "spls", "pam")

# Cross-cell selection frequency for one outcome at one feature level: how often
# each named feature is picked, pooled over every sparse cell in the root.
selection_across_cells <- function(root, level, outcome_id, n = 10) {
  files <- Sys.glob(file.path(
    sweep_root_dir(root), level, "*", "*", "*", "c_data", "results.xlsx"
  ))
  rows <- lapply(files, function(f) {
    sel <- read.xlsx(f, "selection")
    if (!nrow(sel) || !"outcome" %in% names(sel)) {
      return(NULL)
    }
    sel |>
      filter(.data$outcome == outcome_id, .data$model %in% SPARSE_MODELS) |>
      mutate(cell = basename(dirname(dirname(f))))
  })
  bind_rows(rows) |>
    mutate(feature = base_feature(.data$feature)) |>
    group_by(.data$feature) |>
    summarise(
      cells = dplyr::n(), mean_freq = mean(.data$freq), .groups = "drop"
    ) |>
    slice_max(.data$mean_freq, n = n, with_ties = FALSE) |>
    mutate(label = feature_label(.data$feature, level))
}

# Named-driver bars: the interpretable readout the review asks for on every
# composite -- who was selected, how often, under its real name.
ref_named_drivers <- function(root, level, outcome_id, title) {
  d <- selection_across_cells(root, level, outcome_id)
  ggplot(d, aes(.data$mean_freq, stats::reorder(.data$label, .data$mean_freq))) +
    geom_col(fill = SPEC_LEVEL_COLORS[[level]], width = 0.72, alpha = 0.9) +
    geom_text(aes(label = sprintf("%d cells", .data$cells)),
      hjust = -0.15, size = 2.3, colour = "grey25"
    ) +
    scale_x_continuous(
      limits = c(0, 1.15), breaks = seq(0, 1, 0.25), labels = percent
    ) +
    labs(
      x = "mean fold-selection frequency", y = NULL, title = title,
      subtitle = sprintf(
        "Named drivers of %s, pooled across sparse cells",
        outcome_id
      )
    ) +
    FIG_THEME +
    theme(plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE))
}

# The per-cell statistics block the review mandates: metric, raw permutation
# p, how many cells were screened, and the level's leakage status, over the
# null band.
ref_stats_cell <- function(root = "F06_prediction") {
  cells <- root_cells(root)
  top <- cells |> slice_min(.data$perm_p_q2, n = 1, with_ties = FALSE)
  null <- read.xlsx(
    file.path(
      sweep_root_dir(root), top$level, top$config, top$outcome, top$model,
      "c_data", "results.xlsx"
    ),
    "null"
  ) |>
    filter(.data$outcome == top$outcome)

  ggplot(null, aes(.data$q2)) +
    geom_histogram(
      bins = 30, fill = "grey78", colour = "white",
      linewidth = 0.1
    ) +
    geom_vline(xintercept = 0, linetype = "dotted", colour = "grey45") +
    geom_vline(xintercept = top$q2, colour = "#B2182B", linewidth = 1) +
    annotate("label",
      x = -Inf, y = Inf, hjust = -0.05, vjust = 1.15,
      label = sub(
        "  ·  1 of", "\n1 of",
        cell_footer("Q2", top$q2, top$perm_p_q2, root, top$level)
      ),
      size = 2.8, fill = alpha("white", 0.9), lineheight = 0.95
    ) +
    labs(
      x = "permutation-null Q2", y = "count",
      title = sprintf(
        "%s | %s | %s -- %s", top$level, top$config, top$model, top$outcome
      ),
      subtitle = "Observed (red) vs its own null."
    ) +
    FIG_THEME +
    theme(plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE))
}

# Named-feature heatmap: subjects x named drivers, ordered by the phenotype, so
# the signature is legible as biology rather than a feature index.
ref_named_heatmap <- function(root = "F06_prediction", level = "pathways",
                              outcome_id = "d_mcsa", config = "T2") {
  bundle <- pred_load()
  drivers <- selection_across_cells(root, level, outcome_id, n = 10)
  x <- pred_contrast_matrix(
    bundle$feature_sets[[SWEEP_LEVEL_KEY[[level]]]], bundle$meta, config
  )
  y <- pred_outcome(bundle, outcome_id)
  keep <- intersect(rownames(x), names(y)[!is.na(y)])
  x <- x[keep, intersect(drivers$feature, colnames(x)), drop = FALSE]
  z <- scale(x)

  long <- as.data.frame(z) |>
    mutate(subject = rownames(z), pheno = as.numeric(y[keep])) |>
    pivot_longer(-c("subject", "pheno"),
      names_to = "feature",
      values_to = "z"
    ) |>
    mutate(
      label = feature_label(.data$feature, level),
      subject = stats::reorder(.data$subject, .data$pheno)
    )

  ggplot(long, aes(.data$subject, .data$label, fill = .data$z)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "grey96", high = "#B2182B", midpoint = 0,
      name = "z", limits = c(-2.5, 2.5), oob = squish
    ) +
    labs(
      x = sprintf("subjects, ordered by %s", outcome_id), y = NULL,
      title = sprintf("Named drivers -- %s at %s", level, config),
      subtitle = "Rows are real feature names; columns ordered by the phenotype."
    ) +
    FIG_THEME +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
      axis.text.y = element_text(size = 7),
      plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE)
    )
}

# Supplement reference: calibration for the lead cell, plus how selection
# stability decays across the ranked drivers.
ref_calibration <- function(root = "F06_prediction") {
  cells <- root_cells(root)
  top <- cells |> slice_max(.data$q2, n = 1, with_ties = FALSE)
  pr <- read.xlsx(
    file.path(
      sweep_root_dir(root), top$level, top$config, top$outcome, top$model,
      "c_data", "results.xlsx"
    ),
    "predictions"
  ) |>
    filter(.data$outcome == top$outcome) |>
    mutate(group = ifelse(grepl("^HR", .data$subject), "HR", "LR"))
  ggplot(pr, aes(.data$y, .data$pred)) +
    geom_abline(slope = 1, linetype = "dashed", colour = "grey55") +
    geom_smooth(method = "lm", se = TRUE, colour = "grey25", linewidth = 0.5) +
    geom_point(aes(colour = .data$group), size = 2.2) +
    scale_colour_manual(values = GROUP_COLORS, name = NULL) +
    labs(
      x = sprintf("observed %s", top$outcome), y = "predicted",
      title = "Calibration of the lead cell",
      subtitle = sprintf(
        "%s | %s | %s. Dashed = identity; a good forecast tracks it.",
        top$level, top$config, top$model
      )
    ) +
    FIG_THEME +
    theme(
      legend.position = c(0.12, 0.86),
      plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE)
    )
}

ref_stability <- function(root = "F06_prediction", level = "modules",
                          outcome_id = "d_mcsa") {
  d <- selection_across_cells(root, level, outcome_id, n = 13) |>
    mutate(rank = rank(-.data$mean_freq, ties.method = "first"))
  ggplot(d, aes(.data$rank, .data$mean_freq)) +
    geom_hline(yintercept = 0.6, linetype = "dashed", colour = "#B2182B") +
    geom_line(colour = "grey60", linewidth = 0.4) +
    geom_point(size = 2.2, colour = SPEC_LEVEL_COLORS[[level]]) +
    ggrepel::geom_text_repel(aes(label = .data$label),
      size = 2.1,
      max.overlaps = 12, min.segment.length = 0, colour = "grey25"
    ) +
    scale_y_continuous(limits = c(0, 1), labels = percent) +
    labs(
      x = "driver rank", y = "mean selection frequency",
      title = "Selection stability",
      subtitle = paste(
        "Meinshausen-Buhlmann reading: above the dashed line is a candidate,",
        "below it is noise."
      )
    ) +
    FIG_THEME +
    theme(plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE))
}

# QC reference: the upstream that feeds every figure -- module sizes and their
# ORA identity, and the soft-threshold curve behind the module definition.
ref_module_atlas <- function() {
  atlas <- read.xlsx(wgcna_book, "module_atlas")
  lab <- module_ora_label()
  atlas |>
    mutate(
      feature = paste0("ME_", .data$module),
      label = lab$label[match(.data$feature, lab$feature)],
      label = ifelse(is.na(.data$label), .data$module, .data$label)
    ) |>
    ggplot(aes(.data$n_proteins, stats::reorder(.data$label, .data$n_proteins))) +
    geom_col(fill = "#9E9AC8", width = 0.72, alpha = 0.9) +
    geom_text(aes(label = .data$n_proteins),
      hjust = -0.25, size = 2.3,
      colour = "grey25"
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(
      x = "member proteins", y = NULL,
      title = "Module atlas -- what each module is",
      subtitle = "Modules named by their top over-represented term, not by colour."
    ) +
    FIG_THEME +
    theme(plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE))
}

ref_soft_threshold <- function() {
  st <- read.xlsx(wgcna_book, "soft_threshold")
  ggplot(st, aes(.data$Power, .data$signed_r2)) +
    geom_hline(yintercept = 0.8, linetype = "dashed", colour = "#B2182B") +
    geom_line(colour = "grey60") +
    geom_point(size = 2, colour = "#2166AC") +
    labs(
      x = "soft-threshold power", y = "scale-free fit R^2",
      title = "WGCNA soft threshold",
      subtitle = "Dashed = 0.8 scale-free criterion; the chosen power clears it."
    ) +
    FIG_THEME +
    theme(plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE))
}

ref_annotation <- function(title, subtitle) {
  plot_annotation(
    title = title, subtitle = subtitle, tag_levels = "A",
    theme = theme(
      plot.title = element_text(face = "bold", size = FIG_TITLE_SIZE + 3),
      plot.subtitle = element_text(
        face = "italic", size = FIG_SUBTITLE_SIZE, colour = "grey30"
      ),
      plot.tag = element_text(face = "bold", size = FIG_TAG_SIZE)
    )
  )
}

main_ref <-
  (ref_named_drivers("F06_prediction", "proteins", "d_mcsa", "Proteins") |
    ref_named_drivers("F06_prediction", "pathways", "d_mcsa", "Pathways") |
    ref_named_drivers("F06_prediction", "modules", "d_mcsa", "Modules (ORA)")) /
  (ref_named_heatmap() | ref_stats_cell()) +
  ref_annotation(
    "REFERENCE -- main composite idiom (named features + per-cell statistics)",
    paste(
      "A-C: drivers under their real names; modules renamed by ORA, not",
      "colour. D: named-feature heatmap ordered by phenotype.\nE: the",
      "statistics block every panel carries -- metric, raw permutation p,",
      "cells screened, null band, leakage flag. No BH q anywhere."
    )
  )
save_panel(main_ref, file.path(REF_DIR, "REFERENCE_main_idiom"), 330, 200)

supp_ref <- (ref_calibration() | ref_stability()) +
  ref_annotation(
    "REFERENCE -- per-figure validation supplement",
    paste(
      "A: calibration of the lead cell against identity. B: selection",
      "stability across ranked named drivers. Both read at n = 16 with the",
      "honest caveat on the face of the panel."
    )
  )
save_panel(supp_ref, file.path(REF_DIR, "REFERENCE_supp_validation"), 260, 110)

qc_ref <- (ref_module_atlas() | ref_soft_threshold()) +
  ref_annotation(
    "REFERENCE -- shared upstream QC",
    paste(
      "A: module atlas, each module named by its top ORA term and sized by",
      "membership. B: the soft-threshold choice behind those modules. This",
      "figure feeds the module level of every downstream figure."
    )
  )
save_panel(qc_ref, file.path(REF_DIR, "REFERENCE_qc_upstream"), 280, 120)

message("reference panels written to ", REF_DIR)

# Main manuscript composites for the sweep figures (F05 classification, F06
# prediction). Each reads its root roll-up and leaf
# workbooks and assembles one A-F composite in the F02/YvO idiom: a ranked
# grid of per-cell panels (ROC for classification, observed-vs-predicted for
# continuous) carrying each cell's own statistics, the whole-screen
# specification curve, and the feature-level distribution. Descriptive review
# figures; every panel shows the null and never quotes one cell as a result.

pacman::p_load(
  here, dplyr, tidyr, ggplot2, patchwork, forcats, pROC, openxlsx, scales
)
source(here("functions", "shared_style.R"))
source(here("functions", "sweep_grid.R"))
source(here("functions", "sweep_speccurve.R"))

# Leaves carry a phenotype directory between config and model. Classification
# sweeps one outcome and its roll-up has no outcome column, so it resolves to
# the single HR_LR leaf.
leaf_sheet <- function(root, level, config, phenotype, method, sheet) {
  openxlsx::read.xlsx(
    file.path(
      sweep_root_dir(root), level, config, phenotype, method,
      "c_data", "results.xlsx"
    ),
    sheet
  )
}

cell_phenotype <- function(cell) {
  if (is.null(cell$outcome)) "HR_LR" else cell$outcome
}

# Each cell reports at its own best B. Filtering on max(B) over the pooled
# table instead would silently drop every cell swept at a lower B -- the
# documented split runs the fast levels at 0/200/1000 and proteins at 0/200,
# so that filter deletes the entire protein level from the figure with no gap
# and no warning. sweep_rollup.R already solved this; the composites now use
# the same helper.
root_cells <- function(root) {
  openxlsx::read.xlsx(
    file.path(sweep_root_dir(root), "c_data", "results.xlsx"), "all_cells"
  ) |>
    best_b_per_cell()
}

cell_tag <- function(level, config, method) {
  sprintf("%s · %s\n%s", level, config, method)
}

# One filled ROC panel for a single classification cell, recomputed from its
# stored out-of-fold predictions. Border weight and linetype encode nominal
# significance; the fill is the feature-level hue.
roc_panel <- function(root, cell) {
  pr <- leaf_sheet(
    root, cell$level, cell$config, cell_phenotype(cell), cell$model,
    "predictions"
  )
  roc <- pROC::roc(pr$y, pr$pred,
    quiet = TRUE, levels = c(0, 1), direction = "<"
  )
  df <- data.frame(
    fpr = 1 - roc$specificities, tpr = roc$sensitivities
  ) |>
    dplyr::arrange(.data$fpr, .data$tpr)
  sig <- is_lead(cell$estimate, cell$perm_p, "class")
  fill <- SPEC_LEVEL_COLORS[[cell$level]]
  fill_df <- df |>
    group_by(.data$fpr) |>
    summarise(tpr = max(.data$tpr), .groups = "drop")
  ggplot(df, aes(.data$fpr, .data$tpr)) +
    geom_abline(
      slope = 1, linetype = "dashed", colour = "grey60", linewidth = 0.3
    ) +
    geom_area(data = fill_df, fill = fill, alpha = 0.35) +
    geom_step(colour = fill, linewidth = 0.8, direction = "hv") +
    annotate("label",
      x = 0.62, y = 0.22,
      label = sprintf(
        "AUC %.2f\n[%.2f, %.2f]\np = %.3f", cell$estimate, cell$ci_lo,
        cell$ci_hi, cell$perm_p
      ),
      size = 2.2, linewidth = 0, fill = alpha("white", 0.85), lineheight = 0.9
    ) +
    annotate("label",
      x = 0.5, y = 0.05, label = cell_tag(cell$level, cell$config, cell$model),
      size = 2, linewidth = 0, fill = alpha("white", 0.8), lineheight = 0.85,
      fontface = "bold"
    ) +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 1)) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 1)) +
    labs(x = NULL, y = NULL) +
    FIG_THEME +
    theme(
      panel.border = element_rect(
        colour = "black", fill = NA,
        linewidth = if (isTRUE(sig)) 1.4 else 0.3
      ),
      panel.grid = element_blank(),
      axis.text = element_text(size = 6)
    )
}

roc_grid <- function(root, n = 12, ncol = 4) {
  top <- root_cells(root) |>
    arrange(.data$perm_p, dplyr::desc(.data$estimate)) |>
    head(n)
  panels <- lapply(seq_len(nrow(top)), function(i) roc_panel(root, top[i, ]))
  wrap_plots(panels, ncol = ncol)
}

# One observed-vs-predicted panel for a continuous cell, coloured by responder
# group with a fitted line and the Spearman rho.
obs_pred_panel <- function(root, cell) {
  pr <- leaf_sheet(
    root, cell$level, cell$config, cell_phenotype(cell), cell$model,
    "predictions"
  ) |>
    mutate(group = ifelse(grepl("^HR", .data$subject), "HR", "LR"))
  fill <- SPEC_LEVEL_COLORS[[cell$level]]
  sig <- is_lead(cell$q2, cell$perm_p_q2, "cont")
  ggplot(pr, aes(.data$y, .data$pred)) +
    geom_smooth(method = "lm", se = FALSE, colour = fill, linewidth = 0.6) +
    geom_point(aes(colour = group), size = 1.4, alpha = 0.9) +
    scale_colour_manual(values = GROUP_COLORS, guide = "none") +
    annotate("label",
      x = -Inf, y = Inf, hjust = -0.05, vjust = 1.1,
      label = sprintf(
        "Q2 %.2f  p %.3f\n%s · %s · %s", cell$q2, cell$perm_p_q2,
        cell$level, cell$config, cell$model
      ),
      size = 2, linewidth = 0, fill = alpha("white", 0.85), lineheight = 0.9,
      fontface = "bold"
    ) +
    labs(x = sprintf("observed %s", cell$outcome), y = "predicted") +
    FIG_THEME +
    theme(
      panel.border = element_rect(
        colour = "black", fill = NA, linewidth = if (isTRUE(sig)) 1.4 else 0.3
      ),
      axis.title = element_text(size = 6.5), axis.text = element_text(size = 6)
    )
}

obs_pred_grid <- function(root, n = 12, ncol = 4) {
  top <- root_cells(root) |>
    arrange(.data$perm_p_q2, dplyr::desc(.data$q2)) |>
    head(n)
  panels <- lapply(
    seq_len(nrow(top)), function(i) obs_pred_panel(root, top[i, ])
  )
  wrap_plots(panels, ncol = ncol)
}

section_theme <- function() {
  theme(
    plot.title = element_text(face = "bold", size = FIG_TITLE_SIZE),
    plot.subtitle = element_text(
      face = "italic", size = FIG_SUBTITLE_SIZE, colour = "grey30"
    )
  )
}

figure_theme <- function() {
  theme(plot.title = element_text(face = "bold", size = FIG_TITLE_SIZE + 4))
}

build_class_composite <- function(root = "F05_classification") {
  cells <- root_cells(root)
  a <- roc_grid(root, n = 12) +
    plot_annotation(
      title = "A   Top-12 cells by permutation p -- per-cell ROC",
      subtitle = paste(
        "Filled by feature level. Bold border = a lead: beats its null AND",
        "beats chance.\nCells shown are the 12 smallest p, so a thin border",
        "marks one that clears p but not the baseline."
      ),
      theme = section_theme()
    )
  n_lead <- sum(is_lead(cells$estimate, cells$perm_p, "class"), na.rm = TRUE)
  b <- spec_curve_class(
    cells, "B   Specification curve -- every cell vs its null band",
    sprintf(
      paste(
        "%d of %d cells lead; smallest permutation p = %.3f. A feature-free",
        "model under LOSO scores 0, not 0.5, so the 0.5 line is the wrong",
        "baseline for the %d cells at AUC 0."
      ),
      n_lead, nrow(cells), min(cells$perm_p, na.rm = TRUE),
      sum(cells$estimate == 0, na.rm = TRUE)
    )
  )
  wrap_elements(a) / wrap_elements(b) +
    plot_layout(heights = c(1.7, 1.1)) +
    plot_annotation(
      title = "F05 Classification -- can the proteome assign HR vs LR?",
      theme = figure_theme()
    )
}

build_cont_composite <- function(root = "F06_prediction") {
  cells <- root_cells(root)
  a <- obs_pred_grid(root, n = 12) +
    plot_annotation(
      title = "A   Top-12 cells by permutation p -- observed vs predicted",
      subtitle = paste(
        "Points by responder group (HR blue, LR red); line per level.",
        "Bold border = a\nlead: beats its null AND beats predicting the",
        "mean. Cells are the 12 smallest p,\nso a thin border marks one that",
        "clears p on a collapsed null."
      ),
      theme = section_theme()
    )
  b <- spec_curve_cont(
    cells, "B   Specification curve -- every cell vs its null band",
    "Chance = 0. Q^2 clamped at -1. Leads concentrate on d_mcsa."
  )
  wrap_elements(a) / wrap_elements(b) +
    plot_layout(heights = c(1.7, 1.1)) +
    plot_annotation(
      title = "F06 Prediction -- forecasting how much a subject adapts",
      theme = figure_theme()
    )
}

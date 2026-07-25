# Main manuscript composites for the sweep figures (F04 association, F05
# classification, F06 prediction). Each reads its root roll-up and leaf
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

leaf_sheet <- function(root, level, config, method, sheet) {
  openxlsx::read.xlsx(
    file.path(
      sweep_root_dir(root), level, config, method, "c_data", "results.xlsx"
    ),
    sheet
  )
}

root_cells <- function(root) {
  cells <- openxlsx::read.xlsx(
    file.path(sweep_root_dir(root), "c_data", "results.xlsx"), "all_cells"
  )
  dplyr::filter(cells, .data$B == max(.data$B))
}

cell_tag <- function(level, config, method) {
  sprintf("%s | %s\n%s", level, config, method)
}

# One filled ROC panel for a single classification cell, recomputed from its
# stored out-of-fold predictions. Border weight and linetype encode nominal
# significance; the fill is the feature-level hue.
roc_panel <- function(root, cell) {
  pr <- leaf_sheet(root, cell$level, cell$config, cell$model, "predictions")
  roc <- pROC::roc(pr$y, pr$pred, quiet = TRUE, levels = c(0, 1), direction = "<")
  df <- data.frame(
    fpr = 1 - roc$specificities, tpr = roc$sensitivities
  ) |>
    dplyr::arrange(.data$fpr, .data$tpr)
  sig <- cell$perm_p < 0.05
  fill <- SPEC_LEVEL_COLORS[[cell$level]]
  fill_df <- df |>
    group_by(.data$fpr) |>
    summarise(tpr = max(.data$tpr), .groups = "drop")
  ggplot(df, aes(.data$fpr, .data$tpr)) +
    geom_abline(slope = 1, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
    geom_area(data = fill_df, fill = fill, alpha = 0.35) +
    geom_step(colour = fill, linewidth = 0.8, direction = "hv") +
    annotate("label",
      x = 0.62, y = 0.22,
      label = sprintf(
        "AUC %.2f\n[%.2f, %.2f]\np = %.3f", cell$estimate, cell$ci_lo,
        cell$ci_hi, cell$perm_p
      ),
      size = 2.2, label.size = 0, fill = alpha("white", 0.85), lineheight = 0.9
    ) +
    annotate("label",
      x = 0.5, y = 0.05, label = cell_tag(cell$level, cell$config, cell$model),
      size = 2, label.size = 0, fill = alpha("white", 0.8), lineheight = 0.85,
      fontface = "bold"
    ) +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 1)) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 1)) +
    labs(x = NULL, y = NULL) +
    FIG_THEME +
    theme(
      panel.border = element_rect(
        colour = "black", fill = NA,
        linewidth = if (sig) 1.4 else 0.3
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
  pr <- leaf_sheet(root, cell$level, cell$config, cell$model, "predictions") |>
    filter(.data$outcome == cell$outcome) |>
    mutate(group = ifelse(grepl("^HR", .data$subject), "HR", "LR"))
  fill <- SPEC_LEVEL_COLORS[[cell$level]]
  sig <- cell$perm_p_q2 < 0.05
  ggplot(pr, aes(.data$y, .data$pred)) +
    geom_smooth(method = "lm", se = FALSE, colour = fill, linewidth = 0.6) +
    geom_point(aes(colour = group), size = 1.4, alpha = 0.9) +
    scale_colour_manual(values = GROUP_COLORS, guide = "none") +
    annotate("label",
      x = -Inf, y = Inf, hjust = -0.05, vjust = 1.1,
      label = sprintf(
        "Q2 %.2f  p %.3f\n%s | %s | %s", cell$q2, cell$perm_p_q2,
        cell$level, cell$config, cell$model
      ),
      size = 2, label.size = 0, fill = alpha("white", 0.85), lineheight = 0.9,
      fontface = "bold"
    ) +
    labs(x = sprintf("observed %s", cell$outcome), y = "predicted") +
    FIG_THEME +
    theme(
      panel.border = element_rect(
        colour = "black", fill = NA, linewidth = if (sig) 1.4 else 0.3
      ),
      axis.title = element_text(size = 6.5), axis.text = element_text(size = 6)
    )
}

obs_pred_grid <- function(root, n = 12, ncol = 4) {
  top <- root_cells(root) |>
    arrange(.data$perm_p_q2, dplyr::desc(.data$q2)) |>
    head(n)
  panels <- lapply(seq_len(nrow(top)), function(i) obs_pred_panel(root, top[i, ]))
  wrap_plots(panels, ncol = ncol)
}

composite_title <- function(title, subtitle) {
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

build_class_composite <- function(root = "F05_classification") {
  cells <- root_cells(root)
  a <- roc_grid(root, n = 12)
  b <- spec_curve_class(cells, NULL, NULL)$patches$plots[[1]]
  a / b +
    plot_layout(heights = c(1.6, 1)) +
    composite_title(
      "Classification screen -- can the proteome assign HR vs LR?",
      paste(
        "A: top-12 cells by permutation p, filled ROC per feature level, bold",
        "border = nominal p<.05. B: full specification curve. Nested LOSO;",
        "chance = 0.5. One nominal lead at the chance rate -- an honest null."
      )
    )
}

build_cont_composite <- function(root = "F06_prediction") {
  cells <- root_cells(root)
  a <- obs_pred_grid(root, n = 12)
  b <- spec_curve_cont(cells, NULL, NULL)$patches$plots[[1]]
  a / b +
    plot_layout(heights = c(1.6, 1)) +
    composite_title(
      "Continuous-prediction screen -- forecasting how much a subject adapts",
      paste(
        "A: top-12 cells by permutation p, observed vs predicted, points by",
        "responder group. B: full specification curve; Q^2 clamped at -1.",
        "Leads concentrate on d_mcsa but ride on the cohort-imputed proteome."
      )
    )
}

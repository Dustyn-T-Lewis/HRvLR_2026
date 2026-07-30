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

level_factor <- function(x) {
  factor(x, levels = c("proteins", "pathways", "modules"))
}

# Stacked nominal-hit-count bars for the association screen, split by level.
assoc_hit_bars <- function(hits, key, title, subtitle) {
  d <- hits |>
    count(.data[[key]], .data$level, name = "n") |>
    mutate(level = level_factor(.data$level))
  ord <- hits |>
    group_by(.data[[key]]) |>
    summarise(
      tot = dplyr::n(), n_bh = sum(.data$bh < 0.05), .groups = "drop"
    ) |>
    arrange(.data$tot) |>
    mutate(
      lab = ifelse(
        .data$n_bh > 0, sprintf("%d (%d)", .data$tot, .data$n_bh),
        as.character(.data$tot)
      )
    )
  d[[key]] <- factor(d[[key]], levels = ord[[key]])
  ggplot(d, aes(.data$n, .data[[key]], fill = .data$level)) +
    geom_col(width = 0.7, colour = "white", linewidth = 0.2) +
    geom_text(
      data = ord, aes(.data$tot, .data[[key]], label = .data$lab),
      inherit.aes = FALSE, hjust = -0.3, size = 2.5, fontface = "bold",
      colour = "grey20"
    ) +
    scale_fill_manual(values = SPEC_LEVEL_COLORS, name = "level") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.28))) +
    labs(
      x = "associations at p < .05  (parentheses: within-cell BH q < .05)",
      y = NULL, title = title, subtitle = subtitle
    ) +
    FIG_THEME +
    theme(
      legend.position = c(0.78, 0.3),
      plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE)
    )
}

assoc_rate_bars <- function(cell_summary) {
  d <- cell_summary |>
    group_by(config = factor(.data$config, levels = SWEEP_CONFIGS)) |>
    summarise(
      rate = sum(.data$n_nominal) / sum(.data$n_feature), .groups = "drop"
    )
  ggplot(d, aes(.data$config, .data$rate)) +
    geom_col(fill = "grey60", width = 0.7) +
    geom_hline(yintercept = 0.05, linetype = "dashed", colour = "#B2182B") +
    geom_text(aes(label = scales::percent(.data$rate, accuracy = 1)),
      vjust = -0.4, size = 2.3, colour = "grey25"
    ) +
    scale_y_continuous(
      labels = scales::percent, expand = expansion(mult = c(0, 0.15))
    ) +
    labs(
      x = NULL, y = "nominal hit rate",
      title = "Nominal hit rate per config",
      subtitle = "Dashed = 5% chance; excess = structure to test out-of-sample."
    ) +
    FIG_THEME +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),
      plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE)
    )
}

assoc_effect_strip <- function(hits, lead) {
  d <- hits |>
    filter(.data$outcome == lead) |>
    mutate(level = level_factor(.data$level))
  ggplot(d, aes(.data$effect, .data$level, colour = .data$level)) +
    geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
    geom_jitter(height = 0.18, size = 1.2, alpha = 0.7) +
    scale_colour_manual(values = SPEC_LEVEL_COLORS, guide = "none") +
    labs(
      x = sprintf("effect of %s hits", lead), y = NULL,
      title = sprintf("Effect sizes of %s hits", lead),
      subtitle = "Direction of the in-sample associations at nominal p."
    ) +
    FIG_THEME +
    theme(plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE))
}

assoc_euler <- function(hits, lead) {
  sets <- hits |>
    filter(.data$outcome == lead, .data$level == "pathways") |>
    distinct(.data$config, .data$feature) |>
    (\(d) split(d$feature, d$config))()
  sets <- sets[lengths(sets) > 0]
  if (length(sets) < 2) {
    return(
      ggplot() +
        annotate("text", 0, 0,
          label = "too few sets\nto draw",
          size = 3, colour = "grey45", fontface = "italic"
        ) +
        labs(title = sprintf("Pathway overlap -- %s", lead)) +
        theme_void() +
        theme(plot.title = element_text(face = "bold", size = FIG_TITLE_SIZE))
    )
  }
  pal <- scales::alpha(
    c("#E69F00", "#0072B2", "#009E73", "#CC79A7", "#D55E00")[seq_along(sets)],
    0.5
  )
  grob <- grid::grid.grabExpr(print(plot(
    eulerr::euler(sets),
    fills = pal,
    quantities = list(cex = 0.8), labels = list(cex = 0.8)
  )))
  wrap_elements(grob) +
    labs(
      title = sprintf("Pathway overlap -- %s hits", lead),
      subtitle = "Do configs pick the same pathways?"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = FIG_TITLE_SIZE),
      plot.subtitle = element_text(
        face = "italic", size = FIG_SUBTITLE_SIZE, colour = "grey30"
      )
    )
}

build_assoc_composite <- function(root = "F04_association", lead = "d_mcsa") {
  path <- file.path(sweep_root_dir(root), "c_data", "results.xlsx")
  hits <- openxlsx::read.xlsx(path, "hits")
  cs <- openxlsx::read.xlsx(path, "cell_summary")

  a <- assoc_hit_bars(
    hits, "outcome",
    "Nominal hits per outcome", "One adaptation dominates."
  )
  b <- assoc_hit_bars(
    hits, "config",
    "Nominal hits per config", "Mid-training and trajectory lead."
  )
  cc <- assoc_hit_bars(
    hits, "method",
    "Nominal hits per method", "Consistent across the association methods."
  )
  d <- assoc_rate_bars(cs)
  e <- assoc_effect_strip(hits, lead)
  f <- assoc_euler(hits, lead)

  (a | b | cc) / (d | e | f) +
    plot_annotation(
      title = "F04 Association -- in-sample description across the multiverse",
      subtitle = sprintf(
        paste(
          "Every level x config x method x outcome; %d cells. Bars count",
          "nominal p < .05; parentheses count survivors of BH applied within",
          "each cell over that cell's own feature list. A hit is a lead, not a",
          "result: F05/F06 adjudicate out of sample."
        ),
        nrow(cs)
      ),
      tag_levels = "A",
      theme = theme(
        plot.title = element_text(face = "bold", size = FIG_TITLE_SIZE + 4),
        plot.subtitle = element_text(
          face = "italic", size = FIG_SUBTITLE_SIZE, colour = "grey30"
        ),
        plot.tag = element_text(face = "bold", size = FIG_TAG_SIZE)
      )
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

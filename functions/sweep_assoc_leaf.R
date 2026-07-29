# Association leaf driver. One leaf is a feature level x timepoint config x
# method; inside it the method's outcomes are swept and land as one workbook
# sheet each plus a faceted volcano. The runner returns the tidy per-cell tables
# so the root orchestrator can aggregate the whole specification distribution.

pacman::p_load(here, dplyr, tibble, ggplot2, ggrepel, scales)
source(here("functions", "shared_style.R"))
source(here("functions", "sweep_grid.R"))
source(here("functions", "sweep_assoc.R"))
source(here("functions", "pred_features.R"))
source(here("functions", "shared_prediction.R"))

# limma serves both arms, so its leaf carries the six deltas and the HR-vs-LR
# difference; lm and spearman are continuous only, wilcoxon is the group test.
assoc_method_outcomes <- function(method) {
  switch(method,
    limma = c(ADAPT_OUTCOMES, "group_diff"),
    wilcoxon = "group_diff",
    ADAPT_OUTCOMES
  )
}

ASSOC_METHODS <- c("limma", "lm", "spearman", "wilcoxon")

# Compact display name: pathways get their prefix stripped, proteins and modules
# read as-is.
sweep_feature_label <- function(feature, level) {
  if (level == "pathways") clean_pathway_name(feature) else feature
}

run_assoc_leaf <- function(bundle, level, config, method,
                           root = "F04_association") {
  key <- SWEEP_LEVEL_KEY[[level]]
  x0 <- pred_contrast_matrix(bundle$feature_sets[[key]], bundle$meta, config)
  outcomes <- assoc_method_outcomes(method)

  cells <- lapply(outcomes, function(oc) {
    is_group <- oc == "group_diff"
    y_named <- pred_outcome(bundle, if (is_group) "group" else oc)
    al <- align_xy(x0, y_named)
    run_assoc_cell(al$x, al$y, method, is_group) |>
      mutate(outcome = oc, n = length(al$y), .before = 1)
  })
  names(cells) <- outcomes
  tidy <- bind_rows(cells) |>
    mutate(level = level, config = config, method = method, .before = 1)

  dir <- sweep_leaf_dir(root, level, config, method)
  write_sweep_workbook(
    file.path(dir, "c_data", "results.xlsx"),
    c(cells, list(cell_summary = assoc_cell_summary(tidy)))
  )
  panel <- build_assoc_leaf_panel(tidy, level, config, method)
  save_panel(panel, file.path(dir, "b_reports", "panel"),
    width = 210, height = if (length(outcomes) > 1) 140 else 80
  )
  tidy
}

# Per-outcome hit counts for the leaf and the root roll-up.
assoc_cell_summary <- function(tidy) {
  tidy |>
    group_by(.data$level, .data$config, .data$method, .data$outcome) |>
    summarise(
      n_feature = dplyr::n(),
      n_nominal = sum(.data$p < 0.05, na.rm = TRUE),
      n_bh = sum(.data$bh < 0.05, na.rm = TRUE),
      med_abs_effect = stats::median(abs(.data$effect), na.rm = TRUE),
      .groups = "drop"
    )
}

# Faceted effect vs -log10 p volcano, one facet per outcome, nominal hits
# coloured by direction and the strongest few labelled. This is the pre-split
# view: one facet per outcome before split_ fans them into per-phenotype leaves.
build_assoc_leaf_panel <- function(tidy, level, config, method) {
  df <- tidy |>
    filter(is.finite(.data$p), is.finite(.data$effect)) |>
    mutate(
      nlp = -log10(.data$p),
      dir = case_when(
        .data$p < 0.05 & .data$effect > 0 ~ "Up",
        .data$p < 0.05 & .data$effect < 0 ~ "Down",
        TRUE ~ "NS"
      ),
      label = sweep_feature_label(.data$feature, .env$level)
    )
  top <- df |>
    group_by(.data$outcome) |>
    slice_min(.data$p, n = 5, with_ties = FALSE) |>
    ungroup()
  hits <- df |>
    group_by(.data$outcome) |>
    summarise(
      lab = sprintf("p<.05: %d", sum(.data$dir != "NS")), .groups = "drop"
    )

  ggplot(df, aes(.data$effect, .data$nlp)) +
    geom_point(aes(colour = .data$dir), size = 0.8, alpha = 0.6) +
    ggrepel::geom_text_repel(
      data = top, aes(label = .data$label), size = 2, max.overlaps = 8,
      min.segment.length = 0, colour = "grey20"
    ) +
    geom_text(
      data = hits, aes(x = -Inf, y = Inf, label = .data$lab),
      hjust = -0.1, vjust = 1.4, size = 2.3, fontface = "bold",
      colour = "grey30"
    ) +
    facet_wrap(~outcome, scales = "free") +
    scale_colour_manual(values = DIR_COLORS, guide = "none") +
    labs(
      x = "effect", y = expression(-log[10] * " p"),
      title = sprintf("%s · %s association -- %s", level, config, method),
      subtitle = paste(
        "In-sample; nominal p, no correction across the screen.",
        "Proteins/modules cohort-based, optimistic; pathways leakage-free."
      )
    ) +
    FIG_THEME +
    theme(plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE))
}

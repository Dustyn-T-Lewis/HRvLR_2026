# The score heatmap (F5 a, S3 e), after YvO Figure 5A: one row per score (module eigengene or set
# score), a size bar with its name on the left, and 36 columns. The nine limma contrasts are drawn
# as the r-equivalent of t, sign(t) * sqrt(t^2 / (t^2 + df)), so rows compare on one scale. Then
# three relations for fibre area, muscle area and leg press 1RM (T1 level against the T1 trait;
# T1 level against the % gain; T2 - T1 change against the % gain), each as Spearman rho within HR
# and within LR (02_HRvLR/07_Group_Trait_Slopes) and the limma HR - LR slope as an r-equivalent
# (df = people - 4). Solid outline: FDR < alpha_fdr; dashed: p < alpha_nominal only.
# Expects dplyr, tidyr, ggplot2 and patchwork attached and figure_style.R sourced.
r_equivalent <- \(t, df) sign(t) * sqrt(t^2 / (t^2 + df))
heatmap_traits <- c(fcsa_mixed = "fCSA", mcsa = "mCSA", leg_press_1rm = "press")
heatmap_relations <- tibble::tribble(
  ~relation, ~window, ~prefix,
  "Baseline", "baseline_level", "",
  "Predicts gain", "baseline", "pct_",
  "Tracks gain", "training", "pct_"
)
heatmap_contrasts <- tibble::tribble(
  ~contrast, ~section, ~column,
  "Training_HR", "Training", "HR",
  "Training_LR", "Training", "LR",
  "Training_Interaction", "Training", "HR - LR",
  "Acute_HR", "Acute", "HR",
  "Acute_LR", "Acute", "LR",
  "Acute_Interaction", "Acute", "HR - LR",
  "Baseline_HRvLR", "HR - LR", "T1",
  "Trained_HRvLR", "HR - LR", "T2",
  "Acute_HRvLR", "HR - LR", "T3"
)
heatmap_sections <- c("Training", "Acute", "HR - LR", heatmap_relations$relation)

# One row per feature and column, for the features in `contrast_tests` (feature, contrast, t, df,
# p, fdr); `slopes` is group_trait_slopes.rds and `type` its feature_type.
heatmap_cells <- function(contrast_tests, slopes, type) {
  contrasts <- contrast_tests |>
    inner_join(heatmap_contrasts, by = "contrast") |>
    transmute(
      feature, section, column,
      key = contrast, value = r_equivalent(t, df), p, fdr,
      scale = "r-equivalent"
    )
  links <- bind_rows(
    filter(slopes$within_group, feature_type == type) |>
      transmute(window, trait, feature, who = group, value = rho, p, fdr, scale = "Spearman rho"),
    filter(slopes$group_difference, feature_type == type) |>
      transmute(
        window, trait, feature,
        who = "HR - LR", value = r_equivalent(t, n_people - 4), p, fdr,
        scale = "r-equivalent"
      )
  ) |>
    filter(feature %in% contrast_tests$feature)
  phenotype <- heatmap_relations |>
    tidyr::expand_grid(base = names(heatmap_traits)) |>
    mutate(trait = paste0(prefix, base)) |>
    inner_join(links, by = c("window", "trait"), relationship = "many-to-many") |>
    transmute(
      feature,
      section = relation, column = paste(heatmap_traits[base], who),
      key = paste(relation, base, who), value, p, fdr, scale
    )
  cells <- bind_rows(contrasts, phenotype)
  stopifnot(all(count(cells, feature)$n == 36), !anyNA(cells$value))
  cells
}

# `rows`: feature, label and n (the size bar), in display order from the top.
score_heatmap <- function(cells, rows, size_name, title, subtitle) {
  d <- cells |>
    inner_join(rows, by = "feature") |>
    mutate(
      key = factor(key, unique(cells$key)), label = factor(label, rev(rows$label)),
      section = factor(section, heatmap_sections),
      outline = case_when(
        fdr < alpha_fdr ~ sprintf("FDR < %s", alpha_fdr),
        p < alpha_nominal ~ sprintf("p < %s", alpha_nominal),
        .default = NA
      )
    )
  outlines <- rlang::set_names(c("solid", "22"), sprintf(c("FDR < %s", "p < %s"), alpha_fdr))
  grid <- ggplot(d, aes(key, label)) +
    geom_tile(aes(fill = value), filter(d, scale == "r-equivalent"),
      colour = "white",
      linewidth = 0.3
    ) +
    scale_fill_gradient2(
      low = "#762A83", mid = "white", high = "#1B7837", limits = c(-0.5, 0.5),
      oob = scales::squish, name = "r-equivalent\n(contrast, slope)"
    ) +
    ggnewscale::new_scale_fill() +
    geom_tile(aes(fill = value), filter(d, scale == "Spearman rho"),
      colour = "white",
      linewidth = 0.3
    ) +
    scale_fill_gradient2(
      low = "#4393C3", mid = "white", high = "#D6604D", limits = c(-1, 1),
      name = "Spearman rho\n(within group)"
    ) +
    geom_tile(
      aes(linetype = outline), filter(d, !is.na(outline)),
      fill = NA, colour = "black", linewidth = 0.35, width = 0.9, height = 0.9
    ) +
    geom_text(
      aes(label = sprintf("%.2f", value), fontface = if_else(abs(value) >= 0.4, "bold", "plain")),
      size = 1.45
    ) +
    scale_linetype_manual(values = outlines, name = NULL, na.translate = FALSE) +
    scale_x_discrete(labels = tibble::deframe(distinct(d, key, column))) +
    facet_grid(. ~ section, scales = "free_x", space = "free_x") +
    labs(x = NULL, y = NULL, title = title, subtitle = subtitle) +
    figure_theme +
    theme(
      plot.title.position = "panel",
      axis.text.x = element_text(size = 5, angle = 90, hjust = 1, vjust = 0.5),
      axis.text.y = element_blank(), panel.grid = element_blank(),
      strip.text = element_text(face = "bold", size = 5.5), panel.spacing.x = unit(1.2, "mm"),
      legend.position = "bottom", legend.key.width = unit(6, "mm"),
      legend.key.height = unit(2, "mm")
    )
  sizes <- ggplot(rows, aes(n, factor(label, rev(label)))) +
    geom_col(fill = "grey70", width = 0.7) +
    scale_x_reverse(expand = expansion(mult = c(0.05, 0))) +
    labs(x = size_name, y = NULL) +
    figure_theme +
    theme(
      axis.text.y = element_text(size = 5.2), panel.grid.major.y = element_blank(),
      axis.text.x = element_text(size = 4.5)
    )
  # The title sits on the grid; patchwork keeps the same space empty above the size bar.
  (sizes | grid) + plot_layout(widths = c(1, 14))
}

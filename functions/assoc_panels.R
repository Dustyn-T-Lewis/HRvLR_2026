# Panel builders for the F04 association leaves. Every panel states its own
# result on its face: effect with the BH survivor count, so a null reads as
# "small and uncertain" not "absent" (methods doc Part 7). Association is
# descriptive; the panels say so and point to F05 for generalization.

pacman::p_load(ggplot2, dplyr, tidyr, patchwork, ggrepel, scales)

SPACE_LABELS_A <- c(
  proteins = "proteins", pathways = "pathways (singscore)",
  modules = "modules (ME)"
)

# Effect vs -log10 p volcano, BH-significant features coloured by direction and
# the top few labelled by gene. `effect` is a logFC or a slope depending on the
# caller.
assoc_volcano <- function(df, gene_map = NULL, title = NULL,
                          effect_lab = "effect (HR - LR)", n_label = 8) {
  df <- df |>
    filter(is.finite(.data$p), is.finite(.data$effect)) |>
    mutate(
      nlp = -log10(.data$p),
      dir = case_when(
        .data$bh < 0.05 & .data$effect > 0 ~ "Up",
        .data$bh < 0.05 & .data$effect < 0 ~ "Down",
        TRUE ~ "NS"
      ),
      label = if (!is.null(gene_map)) gene_map[.data$feature] else .data$feature
    )
  n_hit <- sum(df$dir != "NS")
  top <- df |>
    slice_min(.data$p, n = n_label, with_ties = FALSE)

  ggplot(df, aes(.data$effect, .data$nlp)) +
    geom_point(aes(colour = dir), size = 1.1, alpha = 0.7) +
    ggrepel::geom_text_repel(
      data = top, aes(label = label), size = 2.3, max.overlaps = 12,
      min.segment.length = 0, colour = "grey20"
    ) +
    scale_colour_manual(values = DIR_COLORS, guide = "none") +
    labs(
      x = effect_lab, y = expression(-log[10] * " p"),
      title = title,
      subtitle = sprintf("BH<.05: %d / %d", n_hit, nrow(df))
    ) +
    FIG_THEME +
    theme(plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE))
}

# One group-DE leaf: three feature-space volcanoes side by side, honest-null
# subtitle across the top.
build_group_leaf_panel <- function(res_by_space, gene_map, contrast, method) {
  vs <- lapply(names(res_by_space), function(sp) {
    assoc_volcano(
      res_by_space[[sp]],
      gene_map = if (sp == "proteins") gene_map else NULL,
      title = SPACE_LABELS_A[[sp]]
    )
  })
  wrap_plots(vs, nrow = 1) +
    plot_annotation(
      title = sprintf("%s HR vs LR — %s", contrast, method),
      subtitle = paste(
        "In-sample differential abundance; BH within contrast.",
        "F05 tests whether any of it generalizes."
      ),
      theme = theme(
        plot.title = element_text(face = "bold", size = FIG_TITLE_SIZE),
        plot.subtitle = element_text(
          face = "italic", size = FIG_SUBTITLE_SIZE, colour = "grey30"
        )
      )
    )
}

# Variance-partition violins across proteins: subject / group / timepoint /
# residual. Subject dominance is the honest headline at n=16.
build_varpart_panel <- function(vp) {
  long <- vp |>
    tidyr::pivot_longer(
      c("subject", "group", "timepoint", "Residuals"),
      names_to = "source", values_to = "frac"
    ) |>
    mutate(source = factor(
      source,
      levels = c("subject", "group", "timepoint", "Residuals")
    ))
  med <- long |>
    group_by(source) |>
    summarise(m = median(frac), .groups = "drop")

  ggplot(long, aes(source, frac, fill = source)) +
    geom_violin(scale = "width", colour = NA, alpha = 0.8) +
    geom_boxplot(width = 0.12, outlier.size = 0.3, fill = "white") +
    geom_text(
      data = med, aes(y = 1.03, label = sprintf("%.0f%%", 100 * m)),
      size = 2.6, fontface = "bold", colour = "grey20"
    ) +
    scale_fill_manual(
      values = c(
        subject = "#2166AC", group = "#B2182B",
        timepoint = "#E69F00", Residuals = "grey70"
      ),
      guide = "none"
    ) +
    scale_y_continuous(limits = c(0, 1.08), labels = scales::percent) +
    labs(
      x = NULL, y = "variance explained",
      title = "Variance partition across proteins",
      subtitle = paste(
        "Each subject's own baseline offset dominates;",
        "group and timepoint are small. Descriptive."
      )
    ) +
    FIG_THEME +
    theme(plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE))
}

# PCA with per-subject trajectories, points coloured by group and shaped by
# timepoint.
build_pca_panel <- function(scores, ve) {
  ggplot(scores, aes(PC1, PC2)) +
    geom_line(aes(group = Subject_ID), colour = "grey80", linewidth = 0.3) +
    geom_point(aes(colour = Group, shape = Timepoint), size = 2) +
    scale_colour_manual(values = GROUP_COLORS) +
    labs(
      x = sprintf("PC1 (%.1f%%)", ve[[1]]),
      y = sprintf("PC2 (%.1f%%)", ve[[2]]),
      title = "Global ordination",
      subtitle = "All samples; subject trajectories in grey."
    ) +
    FIG_THEME +
    theme(plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE))
}

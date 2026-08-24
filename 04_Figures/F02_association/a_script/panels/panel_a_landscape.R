# Panel A: all 60 cells at once. Plotting the smallest adjusted p per cell
# rather than a survivor count, because 59 counts are zero and a bar chart of
# zeros hides how far from the line each cell fell.
if (!exists("summary_tbl")) {
  source(here::here("04_Figures", "F02_association", "a_script", "setup.R"))
}
pacman::p_load(ggplot2, forcats, dplyr)

build_landscape <- function(summary_tbl, tag = "A") {
  d <- summary_tbl |>
    mutate(
      trait = factor(TRAIT_LABELS[.data$phenotype],
        levels = rev(unname(TRAIT_LABELS))
      ),
      level_lab = factor(LEVEL_LABELS[.data$level],
        levels = unname(LEVEL_LABELS)
      ),
      window_lab = WINDOW_LABELS[.data$window],
      neglog = -log10(.data$min_bh)
    )

  p <- ggplot(d, aes(.data$neglog, .data$trait, colour = .data$window_lab)) +
    geom_vline(
      xintercept = -log10(0.05), linetype = "dashed", colour = "grey50",
      linewidth = 0.4
    ) +
    geom_point(size = 2.4, position = position_dodge(0.6)) +
    facet_wrap(~level_lab) +
    scale_colour_manual(values = unname(GROUP_COLORS), name = NULL) +
    labs(
      title = "Every phenotype, every window, every feature level", tag = tag,
      subtitle = paste(
        "Smallest BH-adjusted p per cell; dashed = 0.05.",
        "One of 60 cells clears it"
      ),
      x = expression(-log[10] * "(smallest adjusted p)"), y = NULL
    ) +
    FIG_THEME +
    theme(
      plot.subtitle =
        element_text(size = 6.5, face = "italic", colour = "grey40"),
      legend.position = "top",
      axis.text.y = element_text(size = 6.5),
      panel.grid.major.y = element_blank()
    )

  list(plot = p, audit = d |> dplyr::select(
    level, window, phenotype, n_subjects, n_features, n_nominal, n_bh, min_bh
  ))
}

a_landscape <- build_landscape(summary_tbl)
save_panel(
  a_landscape$plot, file.path(F02_RPT, "panels", "panel_a_landscape"), 190, 105
)
F02_PANELS[["landscape"]] <- a_landscape$plot
F02_AUDIT[["association_summary"]] <- a_landscape$audit

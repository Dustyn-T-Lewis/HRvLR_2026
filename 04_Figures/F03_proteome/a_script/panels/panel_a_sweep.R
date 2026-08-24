# Panel A: how close each label-contrast came. Plotting the smallest adjusted p
# rather than a survivor count, because ten of the twelve cells have a count of
# zero and a bar chart of zeros says nothing about how far from the line they
# fell.
if (!exists("sweep_summary")) {
  source(here::here("04_Figures", "F03_proteome", "a_script", "setup.R"))
}
pacman::p_load(ggplot2, forcats)

build_sweep <- function(sweep_summary, tag = "A") {
  d <- sweep_summary |>
    mutate(
      label_lab = fct_rev(factor(LABEL_NAMES[.data$label],
        levels = unname(LABEL_NAMES)
      )),
      origin = ifelse(.data$internal, "in composite", "independent"),
      neglog = -log10(.data$min_bh)
    )

  p <- ggplot(d, aes(.data$neglog, .data$label_lab,
    colour = .data$contrast,
    shape = .data$origin
  )) +
    geom_vline(
      xintercept = -log10(0.05), linetype = "dashed",
      colour = "grey50", linewidth = 0.4
    ) +
    geom_point(size = 2.8, position = position_dodge(0.5)) +
    scale_colour_manual(values = unname(GROUP_COLORS), name = NULL) +
    scale_shape_manual(values = c(
      `in composite` = 17,
      `independent` = 16
    ), name = NULL) +
    labs(
      title = "How close each label came", tag = tag,
      subtitle = paste(
        "Smallest BH-adjusted p per cell; dashed = 0.05.",
        "Right of the line means at least one survivor"
      ),
      x = expression(-log[10] * "(smallest adjusted p)"), y = NULL
    ) +
    FIG_THEME +
    theme(
      plot.subtitle =
        element_text(size = 6.5, face = "italic", colour = "grey40"),
      legend.position = "top", legend.box = "vertical",
      legend.spacing.y = unit(1, "pt"),
      panel.grid.major.y = element_blank()
    )

  list(plot = p, audit = d |> select(
    label, contrast, internal, n_tested,
    n_nominal, n_bh, min_bh
  ))
}

a_sweep <- build_sweep(sweep_summary)
save_panel(
  a_sweep$plot, file.path(F03_RPT, "panels", "panel_a_sweep"),
  150, 110
)
F03_PANELS[["sweep"]] <- a_sweep$plot
F03_AUDIT[["sweep_summary"]] <- a_sweep$audit

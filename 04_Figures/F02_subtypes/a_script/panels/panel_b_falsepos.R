# Panel B: how often the null data - which has no groups at all - still gets
# assigned more than one component. This is the reason panel A needs a null and
# a bare BIC table would have been misleading.
if (!exists("cluster_cells")) {
  source(here::here("04_Figures", "F02_subtypes", "a_script", "setup.R"))
}
pacman::p_load(ggplot2)

build_falsepos <- function(cluster_cells, tag = "B") {
  d <- cluster_cells |>
    mutate(
      space_lab = SPACE_LABELS[.data$space],
      cell = paste0(.data$n_pc, " PC")
    )

  p <- ggplot(d, aes(
    .data$cell, .data$null_g_above_1,
    fill = .data$space_lab
  )) +
    geom_col(position = position_dodge(0.8), width = 0.7) +
    geom_hline(
      yintercept = 0.05, linetype = "dashed", colour = "grey40",
      linewidth = 0.4
    ) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    scale_fill_manual(values = unname(GROUP_COLORS), name = NULL) +
    labs(
      title = "The method finds clusters in data that has none", tag = tag,
      subtitle = paste(
        "Share of no-cluster null draws assigned more than one component;",
        "dashed = 5%"
      ),
      x = "Principal components retained", y = NULL
    ) +
    FIG_THEME +
    theme(
      plot.subtitle =
        element_text(size = 6.5, face = "italic", colour = "grey40"),
      legend.position = "top"
    )

  list(plot = p, audit = d |> select(space, n_pc, null_g_above_1))
}

b_falsepos <- build_falsepos(cluster_cells)
save_panel(
  b_falsepos$plot, file.path(F02_RPT, "panels", "panel_b_falsepos"),
  160, 90
)
F02_PANELS[["falsepos"]] <- b_falsepos$plot
F02_AUDIT[["null_false_positive_rate"]] <- b_falsepos$audit

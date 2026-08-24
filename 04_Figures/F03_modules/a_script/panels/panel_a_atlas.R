# Panel A: what the twelve modules are. A module result that cannot be named is
# not a result, and three modules carry no enriched term at all, which is worth
# seeing next to the nine that do.
if (!exists("atlas")) {
  source(here::here("04_Figures", "F03_modules", "a_script", "setup.R"))
}
pacman::p_load(ggplot2, forcats, dplyr)

build_atlas <- function(atlas, tag = "A") {
  d <- atlas |>
    mutate(
      named = .data$n_terms > 0,
      label = ifelse(.data$named, module_label(.data$module), .data$module),
      label = fct_reorder(.data$label, .data$n_proteins)
    )

  p <- ggplot(d, aes(.data$n_proteins, .data$label, fill = .data$named)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = .data$n_proteins), hjust = -0.25, size = 2.2) +
    scale_fill_manual(
      values = c(`TRUE` = "#2166AC", `FALSE` = "grey70"),
      labels = c(`TRUE` = "GO-enriched", `FALSE` = "no enriched term"),
      name = NULL
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(
      title = "The twelve co-expression modules", tag = tag,
      subtitle = paste(
        "Proteins per module, named by their strongest GO term.",
        "Universe = 1635 detected proteins"
      ),
      x = "Proteins", y = NULL
    ) +
    FIG_THEME +
    theme(
      plot.subtitle =
        element_text(size = 6.5, face = "italic", colour = "grey40"),
      axis.text.y = element_text(size = 6),
      legend.position = "top",
      panel.grid.major.y = element_blank()
    )

  list(plot = p, audit = atlas)
}

a_atlas <- build_atlas(atlas)
save_panel(
  a_atlas$plot, file.path(F03_RPT, "panels", "panel_a_atlas"), 150, 105
)
F03_PANELS[["atlas"]] <- a_atlas$plot
F03_AUDIT[["module_atlas"]] <- a_atlas$audit

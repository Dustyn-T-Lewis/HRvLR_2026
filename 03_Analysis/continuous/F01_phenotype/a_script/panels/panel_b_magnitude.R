# Panel B: how much every outcome actually changed, across all 16 subjects,
# no group contrast. Mirrors
# categorical/F01_phenotype/a_script/panels/panel_c_forest.R's layout (one
# row per outcome, ordered fibre/muscle/strength) but reports the change
# distribution itself (median, IQR, range) instead of a between-group
# standardized difference -- there is no second group to standardize
# against here. change_magnitude_table() is phenotype_helpers.R's
# group-free counterpart to change_advantage_table().
if (!exists("meta")) {
  source(here::here(
    "03_Analysis", "continuous", "F01_phenotype", "a_script", "setup.R"
  ))
}
pacman::p_load(ggplot2, forcats)

build_magnitude <- function(meta, tag = "B") {
  tbl <- change_magnitude_table(meta)
  top_down <- c(
    "fCSA Type I", "fCSA Type II", "fCSA Mixed", "mCSA",
    "1RM Extension", "1RM Leg Press"
  )
  df <- mutate(tbl, measure = factor(measure, levels = rev(top_down)))

  p <- ggplot(df, aes(median_sd, measure, color = domain)) +
    geom_vline(
      xintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.4
    ) +
    geom_errorbar(
      aes(xmin = min_sd, xmax = max_sd),
      width = 0, linewidth = 0.3, alpha = 0.5
    ) +
    geom_errorbar(
      aes(xmin = q1_sd, xmax = q3_sd),
      width = 0.25, linewidth = 0.6
    ) +
    geom_point(size = 2.6) +
    scale_color_manual(
      values = DOMAIN_COLORS, name = "Domain",
      labels = c(fibre = "Fibre", muscle = "Muscle", strength = "Strength")
    ) +
    labs(
      title = "Change magnitude by domain, all 16 subjects", tag = tag,
      subtitle = paste0(
        "Median (point), IQR (thick bar) and full range (thin bar) of ",
        "T2-T1 change, in SD of that outcome's own change"
      ),
      x = "Post - Pre change (SD of the outcome's own change score)", y = NULL
    ) +
    FIG_THEME +
    theme(
      plot.subtitle =
        element_text(size = 6.5, face = "italic", color = "grey40"),
      legend.position = "inside",
      legend.position.inside = c(0.88, 0.17),
      legend.background =
        element_rect(fill = scales::alpha("white", 0.7), color = NA),
      legend.title = element_text(size = 8, face = "bold"),
      panel.grid.major.y = element_blank()
    )
  list(plot = p, audit = tbl)
}

b_magnitude <- build_magnitude(meta, tag = "B")
save_panel(
  b_magnitude$plot, file.path(F01_RPT, "panels", "panel_b_magnitude"), 150,
  120
)
F01_PANELS[["magnitude"]] <- b_magnitude$plot
F01_AUDIT[["change_magnitude"]] <- b_magnitude$audit

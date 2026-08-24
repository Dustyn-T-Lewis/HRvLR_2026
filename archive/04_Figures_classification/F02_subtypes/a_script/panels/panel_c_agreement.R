# Panel C: the sweep run backwards. Groups drawn from the proteome alone, then
# tested against every phenotype. Nothing was shown the phenotypes before the
# split was made, so a hit here would be a real association rather than a
# restatement - this is the one direction in the project that could name a
# proteome-defined responder group.
#
# Testing the phenotypes on their continuous scale rather than comparing
# binarised labels, because dichotomising the outcome throws away the
# information the test needs.
if (!exists("cluster_phenotype")) {
  source(here::here("04_Figures", "F02_subtypes", "a_script", "setup.R"))
}
pacman::p_load(ggplot2, forcats, dplyr)

# Four views need four hues; the shared palettes top out at three.
VIEW_COLORS <- c("#E69F00", "#0072B2", "#009E73", "#CC79A7")

build_agreement <- function(cluster_phenotype, cluster_calibration, tag = "C") {
  d <- cluster_phenotype |>
    mutate(
      pheno_lab = fct_rev(factor(PHENO_LABELS[.data$phenotype],
        levels = unname(PHENO_LABELS)
      )),
      view_lab = paste0(
        SPACE_LABELS[.data$space], "\n", VIEW_LABELS[.data$view]
      ),
      neglog = -log10(.data$p)
    )

  p <- ggplot(d, aes(.data$neglog, .data$pheno_lab, colour = .data$view_lab)) +
    geom_vline(
      xintercept = -log10(0.05), linetype = "dashed", colour = "grey50",
      linewidth = 0.4
    ) +
    geom_point(size = 2.2, position = position_dodge(0.6), alpha = 0.9) +
    scale_colour_manual(values = VIEW_COLORS, name = NULL) +
    labs(
      title = "Proteome-drawn groups against every phenotype", tag = tag,
      subtitle = paste(
        "Four proteome views x ten phenotypes; dashed = 0.05 unadjusted.",
        "Nothing reaches it"
      ),
      x = expression(-log[10] * "(p)"), y = NULL
    ) +
    FIG_THEME +
    theme(
      plot.subtitle =
        element_text(size = 6.5, face = "italic", colour = "grey40"),
      legend.position = "top", legend.text = element_text(size = 5.5),
      legend.key.size = unit(7, "pt"),
      panel.grid.major.y = element_blank()
    )

  list(plot = p, audit = d |> dplyr::select(
    space, view, phenotype, n, d, p, bh
  ))
}

c_agreement <- build_agreement(cluster_phenotype, cluster_calibration)
save_panel(
  c_agreement$plot, file.path(F02_RPT, "panels", "panel_c_agreement"), 160, 100
)
F02_PANELS[["agreement"]] <- c_agreement$plot
F02_AUDIT[["cluster_vs_phenotype"]] <- c_agreement$audit
F02_AUDIT[["cluster_phenotype_calibration"]] <- cluster_calibration

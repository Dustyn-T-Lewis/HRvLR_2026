# Panel C: why the pathway counts are not reported as findings. fgsea's
# preranked null permutes gene labels, which assumes proteins are exchangeable;
# they are co-regulated here, so a random split of subjects clears BH about as
# often as the real one.
if (!exists("fgsea_null")) {
  source(here::here("04_Figures", "F03_proteome", "a_script", "setup.R"))
}
pacman::p_load(ggplot2)

build_fgsea <- function(fgsea_null, fgsea_calibration, tag = "C") {
  p <- ggplot(fgsea_null, aes(.data$n_hits)) +
    geom_histogram(bins = 24, fill = "grey80", colour = NA) +
    geom_vline(
      xintercept = fgsea_calibration$observed_hits[1],
      colour = "#B2182B", linewidth = 0.8
    ) +
    annotate("text",
      x = fgsea_calibration$observed_hits[1], y = Inf,
      label = sprintf(
        "observed %d\nnull median %d\np = %.2f",
        fgsea_calibration$observed_hits[1],
        fgsea_calibration$null_hits_median[1],
        fgsea_calibration$p_empirical[1]
      ),
      colour = "#B2182B", size = 2.5, hjust = -0.12, vjust = 1.3
    ) +
    labs(
      title = "Pathway significance survives relabelling", tag = tag,
      subtitle = paste(
        "Significant sets under 100 random splits vs the real HR/LR label,",
        "Baseline contrast"
      ),
      x = "Pathway sets at fgsea padj < 0.05", y = "Permutations"
    ) +
    FIG_THEME +
    theme(
      plot.subtitle =
        element_text(size = 6.5, face = "italic", colour = "grey40")
    )

  list(plot = p, audit = fgsea_calibration)
}

c_fgsea <- build_fgsea(fgsea_null, fgsea_calibration)
save_panel(
  c_fgsea$plot, file.path(F03_RPT, "panels", "panel_c_fgsea"),
  150, 90
)
F03_PANELS[["fgsea"]] <- c_fgsea$plot
F03_AUDIT[["fgsea_calibration"]] <- c_fgsea$audit

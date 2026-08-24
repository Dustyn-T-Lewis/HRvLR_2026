# Panel C: every subject placed in the adaptation space, so the sweep's target
# is visible as a continuum rather than as two boxes. The old HR/LR label is
# drawn only to show where a median cut would have fallen through it.
if (!exists("pheno")) {
  source(here::here("04_Figures", "F01_phenotype", "a_script", "setup.R"))
}
pacman::p_load(ggplot2, ggrepel, dplyr)

build_space <- function(pheno, tag = "C") {
  x <- pheno[, CHANGE_TRAITS]
  x <- as.matrix(x[, colSums(!is.na(x)) == nrow(x), drop = FALSE])
  pc <- stats::prcomp(x, center = TRUE, scale. = TRUE)
  ve <- round(100 * pc$sdev^2 / sum(pc$sdev^2))

  d <- tibble::tibble(
    subject = pheno$subject, given = pheno$group_arm,
    PC1 = pc$x[, 1], PC2 = pc$x[, 2]
  )

  p <- ggplot(d, aes(.data$PC1, .data$PC2, colour = .data$given)) +
    geom_hline(yintercept = 0, colour = "grey88", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "grey88", linewidth = 0.3) +
    geom_point(size = 2.6) +
    ggrepel::geom_text_repel(aes(label = .data$subject),
      size = 1.9,
      max.overlaps = 20, seed = 42, show.legend = FALSE
    ) +
    scale_colour_manual(values = GROUP_COLORS, name = "Original label") +
    labs(
      title = "Adaptation is a continuum, not two boxes", tag = tag,
      subtitle = sprintf(
        "PCA of the eight change scores (PC1 %d%%, PC2 %d%%)", ve[1], ve[2]
      ),
      x = sprintf("PC1 (%d%%)", ve[1]), y = sprintf("PC2 (%d%%)", ve[2])
    ) +
    FIG_THEME +
    theme(
      plot.subtitle =
        element_text(size = 6.5, face = "italic", colour = "grey40"),
      legend.position = "inside", legend.position.inside = c(0.86, 0.13),
      legend.background = element_rect(fill = alpha("white", 0.7), colour = NA)
    )

  list(plot = p, audit = d)
}

c_space <- build_space(pheno)
save_panel(
  c_space$plot, file.path(F01_RPT, "panels", "panel_c_space"),
  135, 110
)
F01_PANELS[["space"]] <- c_space$plot
F01_AUDIT[["adaptation_space"]] <- c_space$audit

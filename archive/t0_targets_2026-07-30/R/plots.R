# One theme, one palette, one writer. Every figure is a function returning a ggplot (or a matrix for
# ComplexHeatmap); save_fig() is the only thing that touches disk, sizing from config.

pal_t0 <- function() {
  c(
    HR = "#2166AC", LR = "#B2182B", high = "#2166AC", low = "#B2182B",
    protein = "#4D4D4D", pathway = "#1B7837", module = "#762A83",
    glmnet = "#B2182B", ranger = "#2166AC", null = "grey70"
  )
}

theme_t0 <- function(base_size = 12) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(colour = "grey35"),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "top"
    )
}

save_fig <- function(plot, name, cfg, width = cfg$figures$width, height = cfg$figures$height) {
  dir <- here::here(cfg$paths$outputs, "figures")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(dir, paste0(name, ".", cfg$figures$device))
  if (methods::is(plot, "Heatmap") || methods::is(plot, "HeatmapList")) {
    grDevices::png(path, width = width, height = height, units = "in", res = cfg$figures$dpi)
    ComplexHeatmap::draw(plot)
    grDevices::dev.off()
  } else {
    ggplot2::ggsave(path, plot, width = width, height = height, dpi = cfg$figures$dpi)
  }
  path
}

fig_missingness <- function(qc_miss) {
  ggplot2::ggplot(qc_miss$per_sample, ggplot2::aes(frac, fill = timepoint)) +
    ggplot2::geom_histogram(bins = 20, colour = "white") +
    ggplot2::labs(
      title = "Detection completeness per sample",
      x = "fraction of proteins detected", y = "samples"
    ) +
    theme_t0()
}

fig_cv_performance <- function(results) {
  ggplot2::ggplot(results, ggplot2::aes(roc_auc, wflow_id, colour = model)) +
    ggplot2::geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("p=%.3f", p_perm)),
      hjust = -0.2, size = 3, show.legend = FALSE
    ) +
    ggplot2::facet_grid(target ~ contrast, scales = "free_y") +
    ggplot2::scale_colour_manual(values = pal_t0()) +
    ggplot2::coord_cartesian(xlim = c(0, 1.15)) +
    ggplot2::labs(
      title = "Leave-one-subject-out AUC across the workflow grid",
      subtitle = "Dashed line = chance; p from the label-permutation null",
      x = "LOSO ROC AUC", y = NULL
    ) +
    theme_t0()
}

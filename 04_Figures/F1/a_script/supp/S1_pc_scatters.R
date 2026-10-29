# Every PC of 02_HRvLR/04_Variance against every covariate, one covariate per page, behind S1 d.
# Not a panel: the S1 composite runs it, and it writes S1_pc_scatters.pdf beside the S1 panels.
library(here)
library(dplyr)
library(ggplot2)
source(here("04_Figures", "style", "figure_style.R"))

variance <- "02_HRvLR/04_Variance/c_data/variance.xlsx"
pc_covariates <- read_sheet(variance, "pc_covariates")
pc_scores <- read_sheet(variance, "pc_scores")
covariates <- c(
  timepoint = "time point", group = "HR/LR", subject = "person", blood = "blood index",
  proteins_detected = "proteins detected", weight = "sample weight"
)
components <- unique(pc_covariates$component)
stopifnot(
  setequal(pc_covariates$covariate, names(covariates)), all(components %in% names(pc_scores))
)

scores_long <- tidyr::pivot_longer(
  pc_scores, all_of(components),
  names_to = "component", values_to = "score"
)
scatter_page <- function(covariate) {
  facets <- pc_covariates[pc_covariates$covariate == covariate, ] |>
    mutate(facet = sprintf("%s (%.1f%%), R² %.2f", component, 100 * pc_share, r2))
  d <- scores_long |>
    inner_join(select(facets, component, facet), by = "component") |>
    mutate(x = .data[[covariate]], facet = factor(facet, facets$facet))
  page <- ggplot(d, aes(x, score)) +
    geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey80")
  page <- if (is.numeric(d$x)) {
    page +
      geom_smooth(method = "lm", formula = y ~ x, linewidth = 0.4, colour = "grey30", se = FALSE) +
      geom_point(aes(colour = timepoint, shape = group), size = 1.2)
  } else {
    page +
      geom_boxplot(outliers = FALSE, linewidth = 0.25, width = 0.6) +
      geom_point(
        aes(colour = timepoint, shape = group),
        size = 1.1, position = position_jitter(width = 0.15, seed = 1)
      )
  }
  page +
    scale_colour_manual(values = time_colours, name = NULL) +
    scale_shape_manual(values = group_shapes, name = NULL) +
    facet_wrap(~facet, nrow = 1, scales = "free_y") +
    labs(
      x = covariates[[covariate]], y = "PC score",
      title = sprintf("Principal components against %s", covariates[[covariate]]),
      subtitle = "Blood-corrected values, proteins scaled; R² of the covariate alone on each PC"
    ) +
    figure_theme +
    theme(
      legend.position = "bottom", strip.text = element_text(face = "bold", size = 6),
      axis.text.x = if (covariate == "subject") element_text(angle = 90, hjust = 1, vjust = 0.5)
    )
}
cairo_pdf(
  file.path(figure_dir("S1", "b_reports"), "S1_pc_scatters.pdf"),
  width = 250 / 25.4, height = 80 / 25.4, onefile = TRUE
)
invisible(lapply(lapply(names(covariates), scatter_page), print))
invisible(dev.off())

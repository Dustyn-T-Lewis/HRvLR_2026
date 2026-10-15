library(here)
library(dplyr)
library(ggplot2)
source(here("04_Figures", "style", "figure_style.R"))

inputs <- c(variance = "02_HRvLR/04_Variance/c_data/variance.xlsx")
read_sheet <- \(sheet) readxl::read_excel(here(inputs[["variance"]]), sheet)
pc_covariates <- read_sheet("pc_covariates")
pc_scores <- read_sheet("pc_scores")
holds <- "Panel d: R² of each covariate alone on each of PCs 1 to 5, and each PC's variance share."

covariates <- c(
  timepoint = "time point", group = "HR/LR", subject = "person", blood = "blood index",
  proteins_detected = "proteins detected", weight = "sample weight"
)
components <- unique(pc_covariates$component)
panel_data <- pc_covariates |>
  mutate(
    covariate_label = factor(covariates[covariate], rev(covariates)),
    component_label = sprintf("%s\n%.1f%%", component, 100 * pc_share),
    component_label = factor(component_label, unique(component_label))
  )
stopifnot(
  !anyNA(panel_data$covariate_label), nrow(panel_data) == length(components) * length(covariates),
  all(components %in% names(pc_scores))
)

panel <- ggplot(panel_data, aes(component_label, covariate_label, fill = r2)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(
    aes(label = sprintf("%.2f", r2), colour = r2 > 0.5),
    size = 2, show.legend = FALSE
  ) +
  scale_fill_distiller(
    palette = "Blues", direction = 1, limits = c(0, 1), name = "R²"
  ) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey15")) +
  scale_x_discrete(position = "top") +
  labs(
    x = NULL, y = NULL, title = "What drives the main components",
    subtitle = "R² of each covariate alone; PC share of variance under its name"
  ) +
  figure_theme +
  theme(
    panel.grid.major = element_blank(), legend.position = "bottom",
    legend.title = element_text(vjust = 0.9),
    legend.key.width = unit(6, "mm"), legend.key.height = unit(2, "mm")
  )

ggsave(
  here("04_Figures", "S1", "b_reports", "panels", "d_pc_covariates.pdf"), panel,
  width = 70, height = 70, units = "mm", device = cairo_pdf
)

# Every PC against every covariate, one covariate per page; not part of the composite.
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
  here("04_Figures", "S1", "b_reports", "S1_pc_scatters.pdf"),
  width = 180 / 25.4, height = 75 / 25.4, onefile = TRUE
)
invisible(lapply(lapply(names(covariates), scatter_page), print))
invisible(dev.off())

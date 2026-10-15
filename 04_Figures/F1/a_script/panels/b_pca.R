library(here)
library(dplyr)
library(ggplot2)
library(patchwork)
source(here("04_Figures", "style", "figure_style.R"))

inputs <- c(contrasts = "02_HRvLR/01_Protein_Contrasts/c_data/protein_contrasts.xlsx")
read_sheet <- \(sheet) readxl::read_excel(here(inputs[["contrasts"]]), sheet)
pca <- read_sheet("pca_by_time")
permanova <- read_sheet("permanova")
holds <- "Panel b: PC1 and PC2 per biopsy within each time point, with that time point's PERMANOVA."

group_tests <- permanova |>
  filter(term == "group") |>
  select(timepoint = window, permanova_n = n_samples, permanova_r2 = r2, permanova_p = p)
time_test <- filter(permanova, term == "time")
panel_data <- pca |>
  mutate(group = factor(group, names(group_shapes))) |>
  inner_join(group_tests, by = "timepoint") |>
  select(timepoint, sample_id, subject, group, PC1, PC2, pc1_share, pc2_share, starts_with("perm"))
stopifnot(
  nrow(panel_data) == nrow(pca), nrow(group_tests) == 3, nrow(time_test) == 1,
  all(count(panel_data, timepoint)$n == group_tests$permanova_n)
)

time_names <- c(T1 = "T1 baseline", T2 = "T2 trained", T3 = "T3 acute bout")
time_plot <- function(time) {
  d <- filter(panel_data, timepoint == time)
  axis <- \(i, share) sprintf("PC%d (%.1f%%)", i, 100 * share[1])
  ggplot(d, aes(PC1, PC2)) +
    geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey85") +
    geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey85") +
    geom_point(aes(colour = group, shape = group), size = 1.3) +
    annotate(
      "text",
      x = -Inf, y = Inf, hjust = -0.05, vjust = 1.3, size = 1.9, colour = "grey20",
      label = sprintf(
        "PERMANOVA, group\nR² %.2f, %s", d$permanova_r2[1], label_p(d$permanova_p[1])
      )
    ) +
    scale_colour_manual(values = group_greys, name = NULL, drop = FALSE) +
    scale_shape_manual(values = group_shapes, name = NULL, drop = FALSE) +
    # Headroom above the points for the test.
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.4))) +
    labs(title = time_names[[time]], x = axis(1, d$pc1_share), y = axis(2, d$pc2_share)) +
    figure_theme +
    theme(
      plot.title = element_text(colour = time_colours[[time]]),
      axis.title = element_text(size = 6)
    )
}

# Wrapped as one element with its own title, so the composite tags the panel once.
grid <- wrap_plots(lapply(names(time_names), time_plot)) + guide_area() +
  plot_layout(ncol = 2, guides = "collect")
panel <- wrap_elements(panel = grid) +
  labs(
    title = "Proteome within each time point",
    subtitle = sprintf(
      paste0(
        "PCA of blood-corrected values, proteins scaled\n",
        "PERMANOVA, time across all %d biopsies: R² %.2f, %s"
      ),
      time_test$n_samples, time_test$r2, label_p(time_test$p)
    )
  ) +
  figure_theme

ggsave(
  here("04_Figures", "F1", "b_reports", "panels", "b_pca.pdf"), panel,
  width = 85, height = 80, units = "mm", device = cairo_pdf
)

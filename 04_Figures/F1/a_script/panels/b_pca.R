# F1 b: the proteome within each time point, three PCAs in a row. Each carries the test of HR
# against LR (PERMANOVA) and the two checks for natural clusters (gap statistic k, and how a forced
# two-cluster split matches HR/LR).
library(here)
library(dplyr)
library(ggplot2)
library(patchwork)
source(here("04_Figures", "style", "figure_style.R"))

inputs <- c(
  contrasts = "02_HRvLR/01_Protein_Contrasts/c_data/protein_contrasts.xlsx",
  variance = "02_HRvLR/04_Variance/c_data/variance.xlsx"
)
holds <- paste(
  "Panel b: PC1 and PC2 per biopsy within each time point, with that time point's PERMANOVA,",
  "gap k and forced-split adjusted Rand index."
)
pca <- read_sheet(inputs[["contrasts"]], "pca_by_time")
permanova <- read_sheet(inputs[["contrasts"]], "permanova")
separation <- read_sheet(inputs[["variance"]], "separation")

time_test <- filter(permanova, term == "time")
per_time <- permanova |>
  filter(term == "group") |>
  select(timepoint = window, permanova_n = n_samples, permanova_r2 = r2, permanova_p = p) |>
  inner_join(
    select(separation, timepoint = window, k_gap, ari_group),
    by = "timepoint"
  )
panel_data <- pca |>
  inner_join(per_time, by = "timepoint") |>
  select(
    timepoint, sample_id, subject, group, PC1, PC2, pc1_share, pc2_share,
    permanova_n:ari_group
  )
times <- unique(pca$timepoint)
stopifnot(
  nrow(panel_data) == nrow(pca), nrow(per_time) == length(times), nrow(time_test) == 1,
  all(count(panel_data, timepoint)$n == per_time$permanova_n)
)

time_names <- c(T1 = "T1 baseline", T2 = "T2 trained", T3 = "T3 acute bout")
axis_title <- \(i, share) sprintf("PC%d (%.0f%%)", i, 100 * share)
time_plot <- function(time) {
  d <- filter(panel_data, timepoint == time)
  tests <- with(d[1, ], sprintf(
    "R² %.2f, %s\ngap k = %d; ARI %.2f", permanova_r2, label_p(permanova_p), k_gap, ari_group
  ))
  ggplot(d, aes(PC1, PC2)) +
    geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey85") +
    geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey85") +
    geom_point(aes(colour = group, shape = group), size = 1.5) +
    scale_colour_manual(values = group_greys, name = NULL) +
    scale_shape_manual(values = group_shapes, name = NULL) +
    labs(
      title = time_names[[time]], subtitle = tests,
      x = axis_title(1, d$pc1_share[1]), y = axis_title(2, d$pc2_share[1])
    ) +
    figure_theme +
    theme(
      plot.title = element_text(colour = time_colours[[time]], face = "bold", size = 6.5),
      plot.subtitle = element_text(size = 5.5, colour = "grey25", lineheight = 0.95),
      axis.title = element_text(size = 6)
    )
}

grid <- wrap_plots(lapply(times, time_plot), nrow = 1) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0))
panel <- wrap_elements(panel = grid) +
  labs(
    title = "Proteome within each time point",
    subtitle = with(time_test, sprintf(
      paste0(
        "PCA of blood-corrected values, proteins scaled. Per time point: PERMANOVA for HR/LR;\n",
        "gap k (1 = no natural clusters); ARI, forced two-cluster split vs HR/LR (0 = chance)\n",
        "time across all %d biopsies: PERMANOVA R² %.2f, %s"
      ),
      n_samples, r2, label_p(p)
    ))
  ) +
  figure_theme

save_panel(panel, "F1", "b_pca", 105, 70)

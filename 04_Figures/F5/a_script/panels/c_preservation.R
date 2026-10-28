library(here)
library(dplyr)
library(ggplot2)
source(here("04_Figures", "style", "figure_style.R"))

inputs <- c(preservation = "02_HRvLR/05_Module_Preservation/c_data/preservation.rds")
preservation <- readRDS(here(inputs[["preservation"]]))
holds <- "Panel c: Zsummary of each group's modules tested in the other group."
thresholds <- c(none = 2, strong = 10)
# The axis reaches at least 12, so the strong line stays in view when no module passes it.
z_floor <- 12

panel_data <- preservation$preservation |>
  select(reference, test, module, n_proteins, z_summary, median_rank, full_module, full_share) |>
  arrange(reference, desc(z_summary))

networks <- preservation$group_networks
stopifnot(
  !anyNA(panel_data$z_summary), setequal(panel_data$reference, c("HR", "LR")),
  all(panel_data$reference != panel_data$test),
  all(count(panel_data, reference)$n == networks$n_modules[match(
    count(panel_data, reference)$reference, networks$group
  )])
)

d <- panel_data |>
  mutate(
    facet = factor(sprintf("%s modules in %s", reference, test)),
    row = reorder(paste(reference, module), z_summary),
    label = sprintf("%s → %s (%.0f%%)", module, full_module, 100 * full_share)
  )
panel <- ggplot(d, aes(z_summary, row)) +
  annotate(
    "rect",
    xmin = -Inf, xmax = thresholds[["none"]], ymin = -Inf, ymax = Inf, fill = "grey93"
  ) +
  geom_vline(xintercept = thresholds, linetype = "dashed", linewidth = 0.25, colour = "grey45") +
  geom_point(aes(fill = module, size = n_proteins), shape = 21, colour = "grey25", stroke = 0.25) +
  scale_fill_identity() +
  scale_size_area(max_size = 3, name = "proteins") +
  scale_y_discrete(labels = tibble::deframe(select(d, row, label))) +
  scale_x_continuous(
    limits = c(0, max(z_floor, ceiling(max(d$z_summary)))),
    breaks = sort(c(0, 5, unname(thresholds))),
    expand = expansion(mult = c(0, 0.03))
  ) +
  facet_grid(facet ~ ., scales = "free_y", space = "free_y") +
  labs(
    title = "Module preservation between groups",
    subtitle = sprintf(
      paste(
        "Each group's WGCNA modules tested in the other group's data;\nZsummary < %s none,",
        "> %s strong. Row: group module → full-cohort module\nholding most of its proteins"
      ),
      thresholds[["none"]], thresholds[["strong"]]
    ),
    x = "Zsummary", y = NULL
  ) +
  figure_theme +
  theme(
    axis.text.y = element_text(size = 5.2),
    panel.grid.major.y = element_line(linewidth = 0.15, colour = "grey90"),
    strip.background = element_rect(fill = "grey92", colour = "grey60", linewidth = 0.6),
    strip.text.y = element_text(face = "bold", size = 5.5),
    legend.text = element_text(size = 5.5), legend.title = element_text(size = 6)
  )

save_panel(panel, "F5", "c_preservation", 110, 90)

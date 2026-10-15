library(here)
library(dplyr)
library(ggplot2)
source(here("04_Figures", "style", "figure_style.R"))

inputs <- c(
  design = "04_Figures/assets/design.png",
  contrasts = "02_HRvLR/01_Protein_Contrasts/c_data/protein_contrasts.xlsx"
)
design <- png::readPNG(here(inputs[["design"]]))
biopsies <- readxl::read_excel(here(inputs[["contrasts"]]), "sample_weights")
holds <- "Panel a: the design schematic (assets/design.png) and the people and biopsies it shows."

panel_data <- tibble(
  asset = inputs[["design"]], width_px = dim(design)[2], height_px = dim(design)[1],
  n_people = n_distinct(biopsies$subject),
  n_hr = n_distinct(biopsies$subject[biopsies$group == "HR"]),
  n_lr = n_distinct(biopsies$subject[biopsies$group == "LR"]),
  n_biopsies = nrow(biopsies)
)
stopifnot(
  panel_data$n_hr + panel_data$n_lr == panel_data$n_people,
  setequal(biopsies$timepoint, c("T1", "T2", "T3"))
)

# The schematic is pinned under the title so it sits level with the panel beside it.
panel <- ggplot() +
  annotation_custom(grid::rasterGrob(
    design,
    y = grid::unit(1, "npc"), just = "top", interpolate = TRUE
  )) +
  labs(
    title = "Design and contrasts",
    subtitle = sprintf(
      "%d adults (%d HR, %d LR), %d biopsies; contrast colours carry through b to d",
      panel_data$n_people, panel_data$n_hr, panel_data$n_lr, panel_data$n_biopsies
    )
  ) +
  theme_void(base_size = 7) +
  theme(
    plot.title = figure_theme$plot.title, plot.subtitle = figure_theme$plot.subtitle,
    plot.title.position = "plot"
  )

ggsave(
  here("04_Figures", "F1", "b_reports", "panels", "a_design.pdf"), panel,
  width = 95, height = 62, units = "mm", device = cairo_pdf
)

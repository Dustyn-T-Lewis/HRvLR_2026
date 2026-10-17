# F1 a: the design schematic, with the people and biopsies it shows counted from the saved data.
library(here)
library(dplyr)
library(ggplot2)
source(here("04_Figures", "style", "figure_style.R"))

inputs <- c(
  design = "04_Figures/assets/design.png",
  contrasts = "02_HRvLR/01_Protein_Contrasts/c_data/protein_contrasts.xlsx"
)
holds <- "Panel a: the design schematic (assets/design.png) and the people and biopsies it shows."
design <- png::readPNG(here(inputs[["design"]]))
biopsies <- read_sheet(inputs[["contrasts"]], "sample_weights")

panel_data <- biopsies |>
  summarise(
    n_people = n_distinct(subject), n_hr = n_distinct(subject[group == "HR"]),
    n_lr = n_distinct(subject[group == "LR"]), n_biopsies = n()
  ) |>
  mutate(asset = inputs[["design"]], width_px = ncol(design), height_px = nrow(design))
stopifnot(
  panel_data$n_hr + panel_data$n_lr == panel_data$n_people,
  setequal(biopsies$timepoint, c("T1", "T2", "T3"))
)

# The schematic is pinned under the title so it sits level with the panel beside it.
panel <- ggplot() +
  annotation_custom(grid::rasterGrob(
    design,
    y = unit(1, "npc"), just = "top", interpolate = TRUE
  )) +
  labs(
    title = "Design and contrasts",
    subtitle = with(panel_data, sprintf(
      "%d adults (%d HR, %d LR), %d biopsies; contrast colours carry through the figure",
      n_people, n_hr, n_lr, n_biopsies
    ))
  ) +
  theme_void(base_size = 7) +
  theme(
    plot.title = figure_theme$plot.title, plot.subtitle = figure_theme$plot.subtitle,
    plot.title.position = "plot"
  )

save_panel(panel, "F1", "a_design", 75, 70)

# Figure 1: the design and the whole-proteome picture before any single result. Reads saved
# outputs only. The one computation is the phenotype PCA, a projection with no test; every
# p-value plotted comes from 02_Protein_Contrasts or 03_Set_Contrasts.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})
source(here("00_Setup", "Shared_functions", "analysis_reporting_helpers.R"))
source(here("00_Setup", "Shared_functions", "figure_style.R"))
source(here("00_Setup", "Shared_functions", "figure_panels.R"))
source(here("00_Setup", "Shared_functions", "group_overview_figure.R"))


inputs <- c(
  contrasts = "03_Group_comparisons/HRvLR/02_Protein_Contrasts/c_data/protein_contrasts.xlsx",
  set_tests = "03_Group_comparisons/HRvLR/03_Set_Contrasts/c_data/set_tests.rds",
  gene_sets = "01_Processing/04_Gene_Sets/c_data/gene_sets.rds",
  phenotype = "00_Setup/Input/phenotype.csv",
  design = "05_Figures/HRvLR/F1_Overview/assets/design.png"
)
paths <- map_chr(inputs, here)
if (!all(file.exists(paths))) {
  stop("Missing inputs: ", paste(inputs[!file.exists(paths)], collapse = ", "))
}

# Panel a: the BioRender design panel, pinned under its title so it sits level with b and c.
panel_design <- ggplot() +
  annotation_custom(grid::rasterGrob(
    png::readPNG(paths[["design"]]),
    y = grid::unit(1, "npc"), just = "top", interpolate = TRUE
  )) +
  labs(
    title = "Design and contrasts",
    subtitle = "16 adults, 45 biopsies; contrast colours carry through b to f"
  ) +
  theme_void(base_size = 7) +
  theme(
    plot.title = figure_theme$plot.title, plot.subtitle = figure_theme$plot.subtitle,
    plot.title.position = "plot"
  )

build_group_overview(
  panel_design, inputs, paths,
  groups = c("HR", "LR"), out_dir = here("05_Figures", "HRvLR", "F1_Overview"),
  stem = "f1_overview",
  note = "Panel a is the BioRender design panel (assets/design.png); it reads no data."
)

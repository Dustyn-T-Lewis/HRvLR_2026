# The theme, palettes and panel saver every figure panel sources. Panels set their own titles,
# axes and colours from these; composites set layout, letters and size.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
})

arm_colours <- c(HR = "#2166AC", LR = "#B2182B")
timepoint_colours <- c(T1 = "#E69F00", T2 = "#0072B2", T3 = "#009E73")
level_colours <- c(protein = "#1B9E77", set = "#D95F02", module = "#7570B3")
contrast_order <- c(
  "Training_Interaction", "Acute_Interaction", "Baseline_HRvLR", "Training_HR", "Training_LR",
  "Acute_HR", "Acute_LR", "Trained_HRvLR", "Acute_HRvLR"
)

theme_figure <- function(base_size = 7) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 1),
      plot.subtitle = element_text(colour = "grey30"),
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.key.size = unit(3, "mm"),
      plot.margin = margin(3, 3, 3, 3)
    )
}

# Reads a sheet from a stage 01 to 04 workbook, addressed from the repository root.
read_sheet <- function(path, sheet) readxl::read_excel(here(path), sheet = sheet)

# One PDF and one PNG per panel in the figure's b_reports/panels/, sized in millimetres.
save_panel <- function(plot, figure, name, width = 85, height = 65) {
  dir <- here("05_Figures", figure, "b_reports", "panels")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  for (extension in c("pdf", "png")) {
    ggsave(file.path(dir, paste0(name, ".", extension)), plot,
      width = width, height = height, units = "mm", dpi = 300, bg = "white"
    )
  }
  invisible(plot)
}

# Sources every panel script of a figure, in file order, into its own environment.
load_panels <- function(figure) {
  dir <- here("05_Figures", figure, "a_script", "panels")
  files <- list.files(dir, "^[A-Z]_.*[.]R$")
  set_names(files, sub("[.]R$", "", files)) |>
    map(\(file) source(file.path(dir, file), local = new.env())$value)
}

# Saves the composite as PDF and PNG, and one workbook sheet per panel holding the numbers it
# plots.
save_composite <- function(composite, panels, figure, width, height) {
  out <- here("05_Figures", figure)
  dir.create(file.path(out, "c_data"), recursive = TRUE, showWarnings = FALSE)
  for (extension in c("pdf", "png")) {
    ggsave(file.path(out, "b_reports", paste0(figure, ".", extension)), composite,
      width = width, height = height, units = "mm", dpi = 300, bg = "white"
    )
  }
  writexl::write_xlsx(
    map(panels, "data"),
    file.path(out, "c_data", paste0(figure, "_data.xlsx"))
  )
}

# The chance_expectation sheets of the three levels, one row per level and comparison, with the
# set level's collections pooled.
read_chance <- function() {
  books <- c(
    protein = "02_Differential_Expression/02_Differential/c_data/02_differential.xlsx",
    protein = "02_Differential_Expression/03_Phenotype/c_data/03_phenotype.xlsx",
    set = file.path(
      "03_Pathway_Enrichment/05_classify_and_associate_sets/c_data",
      "05_classify_and_associate_sets.xlsx"
    ),
    module = file.path(
      "04_Network/05_classify_and_associate_modules/c_data",
      "05_classify_and_associate_modules.xlsx"
    )
  )
  imap(books, \(book, level) mutate(read_sheet(book, "chance_expectation"), level = level)) |>
    list_rbind() |>
    summarise(
      across(c(n_tested, n_nominal, n_expected), sum),
      .by = c(level, analysis, comparison)
    ) |>
    mutate(ratio = n_nominal / n_expected, level = factor(level, names(level_colours)))
}

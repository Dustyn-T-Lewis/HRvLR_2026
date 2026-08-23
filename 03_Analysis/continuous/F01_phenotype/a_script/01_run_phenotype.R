# F01 build, continuous tree: render every panel, stitch the composite,
# write one workbook. Mirrors
# categorical/F01_phenotype/a_script/01_run_phenotype.R, two panels instead
# of three and no group-contrast sheets in the workbook.
pacman::p_load(here, patchwork, ggplot2, dplyr, tidyr, tibble, openxlsx)

F01_PANELS <- list()
F01_AUDIT <- list()
source(here(
  "03_Analysis", "continuous", "F01_phenotype", "a_script", "setup.R"
))
source(here("functions", "shared_utils.R"))

clear_dir(F01_RPT)
dir.create(file.path(F01_RPT, "panels"), recursive = TRUE, showWarnings = FALSE)
dir.create(F01_DAT, recursive = TRUE, showWarnings = FALSE)

panel_dir <- here(
  "03_Analysis", "continuous", "F01_phenotype", "a_script", "panels"
)
for (f in c("panel_a_continuum", "panel_b_magnitude")) {
  source(file.path(panel_dir, paste0(f, ".R")))
}

source(here(
  "03_Analysis", "continuous", "F01_phenotype", "a_script", "composite.R"
))

ggsave(file.path(F01_RPT, "F01_phenotype.png"), composite,
  width = 260, height = 180, units = "mm", dpi = 300, bg = "white"
)
ggsave(file.path(F01_RPT, "F01_phenotype.pdf"), composite,
  width = 260, height = 180, units = "mm", device = PDF_DEVICE, bg = "white"
)

composite_summary <- f01_composite_scores(meta) |>
  summarise(
    n = dplyr::n(), mean = mean(value), sd = sd(value),
    median = stats::median(value), min = min(value), max = max(value)
  )

metadata <- tibble::tribble(
  ~field, ~value,
  "figure", "F01 phenotype atlas, continuous tree",
  "design",
  "16 subjects, repeated measures T1/T2 (T3 not used here), no group term",
  "composite axis",
  "every subject's composite hypertrophy score, ranked, no boundary drawn",
  "change magnitude",
  "per-outcome median, IQR and range of T2-T1 change across all subjects",
  "source", "00_input/HRvLR_meta.csv"
)

wb <- createWorkbook()
addWorksheet(wb, "change_magnitude")
writeData(wb, "change_magnitude", F01_AUDIT$change_magnitude)
addWorksheet(wb, "composite_axis")
writeData(wb, "composite_axis", composite_summary, startRow = 1)
writeData(
  wb, "composite_axis", F01_AUDIT$composite_axis,
  startRow = nrow(composite_summary) + 3
)
addWorksheet(wb, "metadata")
writeData(wb, "metadata", metadata)
saveWorkbook(
  wb, file.path(F01_DAT, "F01_phenotype_source_data.xlsx"),
  overwrite = TRUE
)

cat("F01 (continuous) rebuilt: composite, panels, workbook written\n")

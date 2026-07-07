# F01 composite: stitch the rendered phenotype panels into one figure and write
# the source-data workbook. Run after HRvLR_F01_run.R has written the panels.

pacman::p_load(here, patchwork, ggplot2, png, grid, readr, openxlsx)
source(here("04_Figures", "functions", "style.R"))

RPT <- here("04_Figures", "F01", "b_reports")
DAT <- here("04_Figures", "F01", "c_data")

panel <- function(name) {
  path <- file.path(RPT, "panels", paste0(name, ".png"))
  if (!file.exists(path)) stop("Missing panel: ", path, " (run HRvLR_F01_run.R first)")
  ggplot() +
    annotation_custom(rasterGrob(readPNG(path), interpolate = TRUE)) +
    theme_void() +
    theme(plot.margin = margin(2, 2, 2, 2))
}

panels <- c(
  "panel_a_volume_load", "panel_b_1rm_leg", "panel_c_fcsa_mixed",
  "panel_d_fcsa_type1", "panel_e_fcsa_type2", "panel_f_myovision_fcsa",
  "panel_g_hypertrophy_composite", "panel_h_1rm_ext", "panel_i_mcsa"
)

composite <- wrap_plots(lapply(panels, panel), ncol = 3) +
  plot_annotation(
    title = "Phenotype Atlas (HR vs LR)",
    subtitle = "Training volume, strength, and fiber/muscle cross-sectional area",
    caption = "Bars: group mean ± SEM. n = 8 per group; a partial subject drops out where its timepoint is missing.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
      plot.subtitle = element_text(size = 9, color = "grey30", hjust = 0.5),
      plot.caption = element_text(size = 8, color = "grey30", hjust = 0.5)
    )
  )

ggsave(file.path(RPT, "F01_composite.png"), composite,
  width = 470, height = 340, units = "mm", dpi = 300, bg = "white"
)
ggsave(file.path(RPT, "F01_composite.pdf"), composite,
  width = 470, height = 340, units = "mm", device = PDF_DEVICE, bg = "white"
)

# One source-data workbook: an overview sheet, then one sheet per panel's audit table
audit_files <- sort(list.files(DAT, pattern = "^audit_panel_.*\\.csv$", full.names = TRUE))
sheets <- substr(sub("^audit_", "", tools::file_path_sans_ext(basename(audit_files))), 1, 31)
overview <- data.frame(
  sheet = sheets,
  description = gsub("_", " ", sub("^panel_([A-Za-z])_", "Panel \\U\\1: ", sheets, perl = TRUE)),
  stringsAsFactors = FALSE
)
wb <- createWorkbook()
addWorksheet(wb, "overview")
writeData(wb, "overview", overview)
for (i in seq_along(audit_files)) {
  addWorksheet(wb, sheets[i])
  writeData(wb, sheets[i], read_csv(audit_files[i], show_col_types = FALSE))
}
saveWorkbook(wb, file.path(DAT, "F01_source_data.xlsx"), overwrite = TRUE)
cat("F01 composite + workbook saved to", RPT, "\n")

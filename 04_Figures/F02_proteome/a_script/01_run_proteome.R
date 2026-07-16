#!/usr/bin/env Rscript
# F02 build: render every panel (each owns its supplement), stitch the composite,
# write one workbook. The only script in this unit that writes.
pacman::p_load(here, openxlsx)

F02_AUDIT <- list()
source(here("04_Figures", "F02_proteome", "a_script", "setup.R"))
source(here("04_Figures", "shared", "utils.R"))

clear_dir(RPT_DIR)
dir.create(file.path(RPT_DIR, "panels"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(RPT_DIR, "supp"), recursive = TRUE, showWarnings = FALSE)
dir.create(DAT_DIR, recursive = TRUE, showWarnings = FALSE)

panel_dir <- here("04_Figures", "F02_proteome", "a_script", "panels")
for (p in c(
  "panel_a_pca", "panel_b_trajectory", "panel_c_divergence",
  "panel_d_effectsize", "panel_e_dep_counts", "panel_f_pathways"
)) {
  source(file.path(panel_dir, paste0(p, ".R")))
}

# Sheets are ordered by name so the workbook layout is stable across runs.
sheets <- sort(names(F02_AUDIT))
overview <- data.frame(
  sheet = sheets,
  description = gsub("_", " ", sub("^panel_([A-Za-z])_", "Panel \\U\\1: ", sheets, perl = TRUE)),
  stringsAsFactors = FALSE
)

wb <- createWorkbook()
addWorksheet(wb, "overview")
writeData(wb, "overview", overview)
for (s in sheets) {
  addWorksheet(wb, substr(s, 1, 31))
  writeData(wb, substr(s, 1, 31), F02_AUDIT[[s]])
}
saveWorkbook(wb, file.path(DAT_DIR, "F02_proteome_source_data.xlsx"), overwrite = TRUE)

source(here("04_Figures", "F02_proteome", "a_script", "composite.R"))

ggsave(file.path(RPT_DIR, "F02_proteome.png"), composite,
  width = 310, height = 175, units = "mm", dpi = 300, bg = "white"
)
ggsave(file.path(RPT_DIR, "F02_proteome.pdf"), composite,
  width = 235, height = 175, units = "mm", device = PDF_DEVICE, bg = "white"
)

cat("F02 rebuilt: panels, supplements, composite, workbook written\n")

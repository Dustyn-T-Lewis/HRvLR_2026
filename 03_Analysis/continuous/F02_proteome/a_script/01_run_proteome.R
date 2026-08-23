#!/usr/bin/env Rscript
# F02 build, continuous tree: render every panel, stitch the composite,
# write one workbook. Mirrors
# categorical/F02_proteome/a_script/01_run_proteome.R on the reduced
# 3-panel set (see the panel scripts for what was dropped and why).
pacman::p_load(here, openxlsx)

F02_AUDIT <- list()
source(here(
  "03_Analysis", "continuous", "F02_proteome", "a_script", "setup.R"
))
source(here("functions", "shared_utils.R"))

clear_dir(RPT_DIR)
dir.create(file.path(RPT_DIR, "panels"), recursive = TRUE, showWarnings = FALSE)
dir.create(DAT_DIR, recursive = TRUE, showWarnings = FALSE)

panel_dir <- here(
  "03_Analysis", "continuous", "F02_proteome", "a_script", "panels"
)
for (p in c("panel_a_pca", "panel_b_dep_counts", "panel_c_concordance")) {
  source(file.path(panel_dir, paste0(p, ".R")))
}

sheets <- sort(names(F02_AUDIT))
overview <- data.frame(
  sheet = sheets,
  description = gsub(
    "_", " ", sub("^panel_([A-Za-z])_", "Panel \\U\\1: ", sheets, perl = TRUE)
  ),
  stringsAsFactors = FALSE
)

wb <- createWorkbook()
addWorksheet(wb, "overview")
writeData(wb, "overview", overview)
for (s in sheets) {
  addWorksheet(wb, substr(s, 1, 31))
  writeData(wb, substr(s, 1, 31), F02_AUDIT[[s]])
}
saveWorkbook(
  wb, file.path(DAT_DIR, "F02_proteome_source_data.xlsx"),
  overwrite = TRUE
)

source(here(
  "03_Analysis", "continuous", "F02_proteome", "a_script", "composite.R"
))

ggsave(file.path(RPT_DIR, "F02_proteome.png"), composite,
  width = 300, height = 130, units = "mm", dpi = 300, bg = "white"
)
ggsave(file.path(RPT_DIR, "F02_proteome.pdf"), composite,
  width = 300, height = 130, units = "mm", device = PDF_DEVICE, bg = "white"
)

cat("F02 (continuous) rebuilt: panels, composite, workbook written\n")

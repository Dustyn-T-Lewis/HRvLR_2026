#!/usr/bin/env Rscript
# F02 global proteome overview: render the panels into b_reports/panels (+ supp),
# then write the source-data workbook. Panels are stitched into the figure
# manually at one shared height.
pacman::p_load(here, readr, openxlsx)

source(here("04_Figures", "F02", "a_script", "HRvLR_F02_setup.R"))

unlink(setdiff(list.files(RPT_DIR, full.names = TRUE), file.path(RPT_DIR, ".gitkeep")),
  recursive = TRUE
)
dir.create(file.path(RPT_DIR, "panels"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(RPT_DIR, "supp"), recursive = TRUE, showWarnings = FALSE)

for (p in c("A", "B", "C", "F", "H", "I", "J", "M")) {
  source(here("04_Figures", "F02", "a_script", "panels", sprintf("panel_%s.R", p)))
}
source(here("04_Figures", "F02", "a_script", "supp", "panel_G.R"))

audit_files <- sort(list.files(DAT_DIR, pattern = "^audit_panel_.*\\.csv$", full.names = TRUE))
wb <- createWorkbook()
for (f in audit_files) {
  sheet <- substr(sub("^audit_", "", tools::file_path_sans_ext(basename(f))), 1, 31)
  addWorksheet(wb, sheet)
  writeData(wb, sheet, read_csv(f, show_col_types = FALSE))
}
saveWorkbook(wb, file.path(DAT_DIR, "F02_source_data.xlsx"), overwrite = TRUE)

cat("F02 panels + workbook rendered to", RPT_DIR, "\n")

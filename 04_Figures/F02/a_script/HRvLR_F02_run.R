#!/usr/bin/env Rscript
# F02 global proteome overview: render the panels into b_reports/panels (+ supp),
# then write the source-data workbook. Panels are stitched into the figure
# manually at one shared height.
pacman::p_load(here, readr, openxlsx)

source(here("04_Figures", "F02", "a_script", "HRvLR_F02_setup.R"))

unlink(setdiff(list.files(RPT_DIR, full.names = TRUE), file.path(RPT_DIR, ".gitkeep")),
  recursive = TRUE
)
unlink(list.files(DAT_DIR, pattern = "^audit_panel_.*\\.csv$", full.names = TRUE))
dir.create(file.path(RPT_DIR, "panels"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(RPT_DIR, "supp"), recursive = TRUE, showWarnings = FALSE)

for (p in c("A", "B", "C", "D", "E", "F")) {
  source(here("04_Figures", "F02", "a_script", "panels", sprintf("panel_%s.R", p)))
}
for (s in c("direction", "C", "M", "CAP", "G", "I", "J")) {
  source(here("04_Figures", "F02", "a_script", "supp", sprintf("panel_%s.R", s)))
}

audit_files <- sort(list.files(DAT_DIR, pattern = "^audit_panel_.*\\.csv$", full.names = TRUE))
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
saveWorkbook(wb, file.path(DAT_DIR, "F02_source_data.xlsx"), overwrite = TRUE)

cat("F02 panels + workbook rendered to", RPT_DIR, "\n")

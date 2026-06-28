# F01 phenotype atlas: render panels (a-i) into b_reports/panels, then stitch.
pacman::p_load(here)

RPT <- here("04_Figures", "F01", "b_reports")
unlink(setdiff(list.files(RPT, full.names = TRUE), file.path(RPT, ".gitkeep")),
  recursive = TRUE
)
dir.create(file.path(RPT, "panels"), recursive = TRUE, showWarnings = FALSE)

for (p in LETTERS[1:9]) {
  source(here("04_Figures", "F01", "a_script", "panels", sprintf("panel_%s.R", p)))
}

source(here("04_Figures", "F01", "a_script", "HRvLR_F01_composite.R"))

cat("F01 phenotype panels + composite rendered to", RPT, "\n")

# F02 build. The only script in this unit that writes.
pacman::p_load(here, ggplot2, patchwork, openxlsx, dplyr, tibble)

F02_PANELS <- list()
F02_AUDIT <- list()

a_script <- here("04_Figures", "F02_association", "a_script")
source(file.path(a_script, "setup.R"))
source(here("functions", "shared_utils.R"))

clear_dir(F02_RPT)
dir.create(file.path(F02_RPT, "panels"), recursive = TRUE, showWarnings = FALSE)
dir.create(F02_DAT, recursive = TRUE, showWarnings = FALSE)

for (p in c("panel_a_landscape", "panel_b_hit", "panel_c_robust")) {
  source(file.path(a_script, "panels", paste0(p, ".R")))
}

if (!is.null(sweep_calibration)) {
  F02_AUDIT[["sweep_calibration"]] <- sweep_calibration
  F02_AUDIT[["cell_calibration"]] <- confirmation
}

metadata <- tibble::tribble(
  ~field, ~value,
  "figure", "F02 continuous proteome-to-phenotype association",
  "design", paste(
    "three feature levels by two windows by ten phenotypes, 60 cells.",
    "Each subject contributes one column: the feature change over the",
    "window, regressed on that subject's adaptation"
  ),
  "no baseline", paste(
    "a baseline association compares levels between people, a different",
    "question from whether a change tracks a change; V1 already tested the",
    "baseline form across 54 cells without promoting anything"
  ),
  "no blocking", paste(
    "one row per subject means no repeated measures inside the fit, so no",
    "duplicateCorrelation: the within-subject structure is spent forming",
    "the difference"
  ),
  "estimator", paste(
    "limma with a continuous predictor, coef 2. The moderated variance is",
    "the reason to prefer it to a per-feature lm at n = 14"
  ),
  "multiplicity", paste(
    "BH within each cell, never across the 60. The phenotypes are not ten",
    "independent questions: three fibre-area measures share r > 0.9 and the",
    "fibre counts run inverse to them, leaving about five independent axes"
  ),
  "robustness", paste(
    "any cell clearing BH is refit dropping each subject in turn and",
    "checked against a rank correlation, because a fit at n = 15 can clear",
    "a threshold on two or three points"
  ),
  "source", "03_Features/c_data"
)

sheets <- sort(names(F02_AUDIT))
wb <- createWorkbook()
addWorksheet(wb, "overview")
writeData(wb, "overview", data.frame(sheet = sheets))
for (s in sheets) {
  addWorksheet(wb, s)
  writeData(wb, s, F02_AUDIT[[s]])
}
addWorksheet(wb, "metadata")
writeData(wb, "metadata", metadata)
saveWorkbook(wb, file.path(F02_DAT, "F02_association_source_data.xlsx"),
  overwrite = TRUE
)

source(file.path(a_script, "composite.R"))

ggsave(file.path(F02_RPT, "F02_association.png"), composite,
  width = 300, height = 210, units = "mm", dpi = 300, bg = "white"
)
ggsave(file.path(F02_RPT, "F02_association.pdf"), composite,
  width = 300, height = 210, units = "mm", device = PDF_DEVICE, bg = "white"
)

message("F02 rebuilt")

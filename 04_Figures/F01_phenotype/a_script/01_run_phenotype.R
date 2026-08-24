# F01 build. The only script in this unit that writes.
pacman::p_load(here, ggplot2, patchwork, openxlsx, dplyr, tibble)

F01_PANELS <- list()
F01_AUDIT <- list()

a_script <- here("04_Figures", "F01_phenotype", "a_script")
source(file.path(a_script, "setup.R"))
source(here("functions", "shared_utils.R"))

clear_dir(F01_RPT)
dir.create(file.path(F01_RPT, "panels"),
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create(F01_DAT, recursive = TRUE, showWarnings = FALSE)

for (p in c("panel_a_change", "panel_b_continuum", "panel_c_separation")) {
  source(file.path(a_script, "panels", paste0(p, ".R")))
}

F01_AUDIT[["composite_structure"]] <- composite_structure
F01_AUDIT[["composite_modality"]] <- composite_modality

metadata <- tibble::tribble(
  ~field, ~value,
  "figure", "F01 responder phenotype and label construction",
  "n", "16 subjects, complete on every outcome except 1RM leg extension (15)",
  "effect size", paste(
    "mean or group difference divided by the SD of the change score;",
    "95% CI from the same t interval"
  ),
  "composite", paste(
    "comp_hypertrophy arrives from the source spreadsheet without a",
    "stated formula; the five outcomes reconstruct it at r2 = 0.977"
  ),
  "label", paste(
    "HR/LR is the exact median cut of comp_hypertrophy: the top eight",
    "are HR without exception"
  ),
  "bimodality", paste(
    "mclust BIC selects two components; bootstrap LRT over 999",
    "replicates gives p = 0.035, and that solution recovers the given",
    "label exactly (ARI = 1)"
  ),
  "internal flag", paste(
    "an outcome with r2 >= 0.20 against the composite is an ingredient",
    "of it and cannot corroborate the label cut from it"
  ),
  "not claimed", paste(
    "this figure reports the construction and the separations; it draws",
    "no inference about the proteomic result"
  ),
  "source", "03_Features/01_Responsiveness/c_data/01_label_audit.xlsx"
)

sheets <- sort(names(F01_AUDIT))
wb <- createWorkbook()
addWorksheet(wb, "overview")
writeData(wb, "overview", data.frame(sheet = sheets))
for (s in sheets) {
  addWorksheet(wb, s)
  writeData(wb, s, F01_AUDIT[[s]])
}
addWorksheet(wb, "metadata")
writeData(wb, "metadata", metadata)
saveWorkbook(wb, file.path(F01_DAT, "F01_phenotype_source_data.xlsx"),
  overwrite = TRUE
)

source(file.path(a_script, "composite.R"))

ggsave(file.path(F01_RPT, "F01_phenotype.png"), composite,
  width = 290, height = 190, units = "mm", dpi = 300, bg = "white"
)
ggsave(file.path(F01_RPT, "F01_phenotype.pdf"), composite,
  width = 290, height = 190, units = "mm", device = PDF_DEVICE, bg = "white"
)

message("F01 rebuilt")

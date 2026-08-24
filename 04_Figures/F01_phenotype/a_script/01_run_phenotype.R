# F01 build. The only script in this unit that writes.
pacman::p_load(here, ggplot2, patchwork, openxlsx, dplyr, tibble)

F01_PANELS <- list()
F01_AUDIT <- list()

a_script <- here("04_Figures", "F01_phenotype", "a_script")
source(file.path(a_script, "setup.R"))
source(here("functions", "shared_utils.R"))

clear_dir(F01_RPT)
dir.create(file.path(F01_RPT, "panels"), recursive = TRUE, showWarnings = FALSE)
dir.create(F01_DAT, recursive = TRUE, showWarnings = FALSE)

for (p in c("panel_a_change", "panel_b_structure", "panel_c_space")) {
  source(file.path(a_script, "panels", paste0(p, ".R")))
}

metadata <- tibble::tribble(
  ~field, ~value,
  "figure", "F01 the adaptation phenotypes, as continuous outcomes",
  "n", paste(
    "16 subjects; 1RM leg extension has one missing value and drops that",
    "subject from its own cell only"
  ),
  "effect size", paste(
    "mean change divided by the SD of the change score, with the 95% CI",
    "from the same t interval"
  ),
  "volume_load", paste(
    "a total, not a change, so it carries no row in panel A; it appears in",
    "the correlation structure and is tested like the rest"
  ),
  "fibre counts", paste(
    "the MyoVision columns are fibre counts, not areas, despite the fCSA in",
    "their meta names; they move against area because larger fibres pack",
    "fewer into the imaged field"
  ),
  "no group split", paste(
    "the original HR/LR label appears in panel C only to show where a median",
    "cut would have fallen; nothing in this project now conditions on it"
  ),
  "source", "00_input/c_data/phenotype.csv"
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
  width = 290, height = 200, units = "mm", dpi = 300, bg = "white"
)
ggsave(file.path(F01_RPT, "F01_phenotype.pdf"), composite,
  width = 290, height = 200, units = "mm", device = PDF_DEVICE, bg = "white"
)

message("F01 rebuilt")

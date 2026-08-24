# F03 build. The only script in this unit that writes.
pacman::p_load(here, ggplot2, patchwork, openxlsx, dplyr, tibble)

F03_PANELS <- list()
F03_AUDIT <- list()

a_script <- here("04_Figures", "F03_modules", "a_script")
source(file.path(a_script, "setup.R"))
source(here("functions", "shared_utils.R"))

clear_dir(F03_RPT)
dir.create(file.path(F03_RPT, "panels"), recursive = TRUE, showWarnings = FALSE)
dir.create(F03_DAT, recursive = TRUE, showWarnings = FALSE)

for (p in c("panel_a_atlas", "panel_b_heatmap", "panel_c_calibration")) {
  source(file.path(a_script, "panels", paste0(p, ".R")))
}

F03_AUDIT[["hit_robustness"]] <- robustness

metadata <- tibble::tribble(
  ~field, ~value,
  "figure", "F03 module atlas and module-phenotype association",
  "modules", paste(
    "12 WGCNA modules plus grey, defined on within-subject-centred abundance",
    "and scored on raw, so subject identity cannot drive them while",
    "between-subject contrasts survive"
  ),
  "annotation", paste(
    "clusterProfiler::enrichGO over BP, CC and MF with the universe set to",
    "the 1635 detected proteins, not the genome. Labelling only: nothing",
    "downstream gates on an enrichment q"
  ),
  "grid", paste(
    "12 modules x 10 phenotypes x 6 windows = 720 tests, BH within each of",
    "the 60 module cells"
  ),
  "result", paste(
    "three module-phenotype cells clear BH inside themselves, all against",
    "change in whole-muscle CSA. None survives BH across the 180-cell sweep,",
    "and the sweep produced 3 such cells where 9 are expected by chance"
  ),
  "greenyellow", paste(
    "the extracellular-matrix module, and the only cell to pass rank",
    "correlation, leave-one-subject-out and adjustment for biopsy",
    "composition. It still fails the sweep-level correction, and its fit",
    "rests on LR_S14, the extreme subject on both axes"
  ),
  "composition", paste(
    "at T2 the myofibre fraction correlates -0.81 with change in",
    "whole-muscle CSA and blood +0.71, so level-window associations with",
    "that phenotype are partly about biopsy content"
  ),
  "source", "03_Features/c_data"
)

sheets <- sort(names(F03_AUDIT))
wb <- createWorkbook()
addWorksheet(wb, "overview")
writeData(wb, "overview", data.frame(sheet = sheets))
for (s in sheets) {
  addWorksheet(wb, s)
  writeData(wb, s, F03_AUDIT[[s]])
}
addWorksheet(wb, "metadata")
writeData(wb, "metadata", metadata)
saveWorkbook(wb, file.path(F03_DAT, "F03_modules_source_data.xlsx"),
  overwrite = TRUE
)

source(file.path(a_script, "composite.R"))

ggsave(file.path(F03_RPT, "F03_modules.png"), composite,
  width = 330, height = 200, units = "mm", dpi = 300, bg = "white"
)
ggsave(file.path(F03_RPT, "F03_modules.pdf"), composite,
  width = 330, height = 200, units = "mm", device = PDF_DEVICE, bg = "white"
)

message("F03 rebuilt")

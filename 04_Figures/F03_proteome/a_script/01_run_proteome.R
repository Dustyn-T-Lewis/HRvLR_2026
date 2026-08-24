# F03 build. The only script in this unit that writes.
pacman::p_load(here, ggplot2, patchwork, openxlsx, dplyr, tibble)

F03_PANELS <- list()
F03_AUDIT <- list()

a_script <- here("04_Figures", "F03_proteome", "a_script")
source(file.path(a_script, "setup.R"))
source(here("functions", "shared_utils.R"))

clear_dir(F03_RPT)
dir.create(file.path(F03_RPT, "panels"),
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create(F03_DAT, recursive = TRUE, showWarnings = FALSE)

for (p in c("panel_a_sweep", "panel_b_confirm", "panel_c_fgsea")) {
  source(file.path(a_script, "panels", paste0(p, ".R")))
}

metadata <- tibble::tribble(
  ~field, ~value,
  "figure", "F03 label sweep and its calibration",
  "design", paste(
    "ten candidate labels resolving to eight distinct splits, by three",
    "contrasts. Baseline asks whether the groups differed before training;",
    "the two interactions whether they diverged over training and over the",
    "acute bout"
  ),
  "estimator", paste(
    "limma on the six label-by-timepoint cells, subject as the blocking",
    "factor with a consensus duplicateCorrelation; identical to V1's fit to",
    "5e-15 under the given label"
  ),
  "multiplicity", paste(
    "BH within each label and contrast, never across the twelve cells:",
    "five of six labels are splits of correlated outcomes on the same 16",
    "subjects"
  ),
  "no blood covariate", paste(
    "V1 adjusted the acute contrasts for the blood index because T3 biopsies",
    "are bloodier; neither contrast here touches T3"
  ),
  "result", paste(
    "three protein-contrast survivors across 24 cells, all at Baseline.",
    "Two cells cleared BH against 2.6 expected by chance (p = 0.75), and",
    "neither survives a 999-permutation subject-label null"
  ),
  "pathways", paste(
    "fgsea returned 1398 set-contrast hits, but a random split of the same",
    "subjects returns as many, so its padj carries no information here;",
    "singscore through the protein estimator returned zero"
  ),
  "source", paste(
    "03_Features/04_Proteins/c_data and 03_Features/05_Pathways/c_data"
  )
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
saveWorkbook(wb, file.path(F03_DAT, "F03_proteome_source_data.xlsx"),
  overwrite = TRUE
)

source(file.path(a_script, "composite.R"))

ggsave(file.path(F03_RPT, "F03_proteome.png"), composite,
  width = 300, height = 195, units = "mm", dpi = 300, bg = "white"
)
ggsave(file.path(F03_RPT, "F03_proteome.pdf"), composite,
  width = 300, height = 195, units = "mm", device = PDF_DEVICE, bg = "white"
)

message("F03 rebuilt")

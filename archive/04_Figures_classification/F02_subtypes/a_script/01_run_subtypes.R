# F02 build. The only script in this unit that writes.
pacman::p_load(here, ggplot2, patchwork, openxlsx, dplyr, tibble)

F02_PANELS <- list()
F02_AUDIT <- list()

a_script <- here("04_Figures", "F02_subtypes", "a_script")
source(file.path(a_script, "setup.R"))
source(here("functions", "shared_utils.R"))

clear_dir(F02_RPT)
dir.create(file.path(F02_RPT, "panels"),
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create(F02_DAT, recursive = TRUE, showWarnings = FALSE)

for (p in c("panel_a_null", "panel_b_falsepos", "panel_c_agreement")) {
  source(file.path(a_script, "panels", paste0(p, ".R")))
}

metadata <- tibble::tribble(
  ~field, ~value,
  "figure", "F02 unsupervised subtype discovery",
  "question", paste(
    "whether the baseline proteome carries its own two-group structure,",
    "fitted without reference to any label"
  ),
  "engine", paste(
    "mclust::Mclust over G = 1:4, BIC-selected, on principal components of",
    "the feature space; G = 1 is a selectable outcome"
  ),
  "null", paste(
    "999 draws from a single multivariate Gaussian with the observed mean",
    "and covariance: same n, same dimension, same correlations, no clusters"
  ),
  "why a null", paste(
    "at 15 subjects a mixture model assigns more than one component to",
    "structureless data in 47 to 77 percent of draws"
  ),
  "dimension", paste(
    "swept over 2, 3 and 4 components and all reported; choosing the one",
    "that gave the best answer is the failure this sweep prevents"
  ),
  "gate", paste(
    "shut: no cell beat its null at p < 0.05, so no seventh candidate label",
    "was written and the stage-04 sweep stayed at six"
  ),
  "source", "03_Features/03_Subtypes/c_data/01_subtypes.xlsx"
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
saveWorkbook(wb, file.path(F02_DAT, "F02_subtypes_source_data.xlsx"),
  overwrite = TRUE
)

source(file.path(a_script, "composite.R"))

ggsave(file.path(F02_RPT, "F02_subtypes.png"), composite,
  width = 300, height = 200, units = "mm", dpi = 300, bg = "white"
)
ggsave(file.path(F02_RPT, "F02_subtypes.pdf"), composite,
  width = 300, height = 200, units = "mm", device = PDF_DEVICE, bg = "white"
)

message("F02 rebuilt")

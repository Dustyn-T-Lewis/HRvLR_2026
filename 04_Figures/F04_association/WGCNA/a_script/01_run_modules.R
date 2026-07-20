#!/usr/bin/env Rscript
# F04 build: the module atlas, what the modules are, and the two null read-outs.
# The only script in this unit that writes.
pacman::p_load(here, patchwork, ggplot2, openxlsx, dplyr, tibble)

F04_PANELS <- list()
F04_AUDIT <- list()

a_script <- here("04_Figures", "F04_association", "WGCNA", "a_script")
source(file.path(a_script, "setup.R"))
source(here("functions", "shared_utils.R"))

clear_dir(RPT_DIR)
dir.create(file.path(RPT_DIR, "panels"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(RPT_DIR, "supp"), recursive = TRUE, showWarnings = FALSE)
dir.create(DAT_DIR, recursive = TRUE, showWarnings = FALSE)

for (p in c(
  "panel_a_atlas", "panel_b_pathways", "panel_c_fry", "panel_d_phenotype"
)) {
  source(file.path(a_script, "panels", paste0(p, ".R")))
}

source(file.path(a_script, "composite.R"))

ggsave(file.path(RPT_DIR, "WGCNA.png"), composite,
  width = 300, height = 340, units = "mm", dpi = 200, bg = "white"
)
ggsave(file.path(RPT_DIR, "WGCNA.pdf"), composite,
  width = 300, height = 340, units = "mm", device = PDF_DEVICE, bg = "white"
)

F04_AUDIT[["module_membership"]] <- tibble::enframe(
  mods,
  name = "gene", value = "module"
)
F04_AUDIT[["module_eigengene"]] <- me_long

# Eigengene matrix in long form for the F06 prediction feature layer: bare module id
# (moduleEigengenes prefixes "ME"), the sample, and the score. F06 reads group_id/sample_id/ME.
readr::write_csv(
  me_long |> transmute(group_id = sub("^ME", "", module), sample_id, ME),
  file.path(DAT_DIR, "wgcna_eigengene.csv")
)

metadata <- tibble::tribble(
  ~field, ~value,
  "figure", "F04 co-expression modules (HRvLR)",
  "engine", sprintf(
    "WGCNA blockwiseModules: signed, bicor (maxPOutliers = 0.05), power = %d (signed R2 = %.3f, mean k = %.1f), deepSplit 2, minModuleSize 30, mergeCutHeight 0.15, pamRespectsDendro TRUE",
    wg$power, wg$r2, wg$mean_k
  ),
  "module definition", "within-subject-centred abundance (removes subject identity as a driver)",
  "module scoring", "eigengenes computed on RAW abundance, so the between-subject HR/LR contrast survives to be tested",
  "why", "built on raw abundance the modules encode subject identity (ICC up to 0.51); that artifact produced the earlier baseline-prediction q2 of 0.713",
  "module movement", "limma::fry, modules as the index, on the 03_DEP design with subject block and duplicateCorrelation",
  "not fgsea", "the previous NES panel ranked the modules by moderated t from the matrix that defined them - circular, and it returned padj = 7e-33",
  "result", "modules are biologically coherent; none moves with responder status and none predicts adaptation",
  "source", "02_Normalization/imputation (missForest), 00_input/c_data/phenotype.csv, 03_DEP/a_non_imputed"
)

sheets <- c(
  "module_membership", "module_atlas", "module_ora", "module_fry",
  "module_prediction", "module_phenotype_lmm", "module_eigengene", "soft_threshold"
)
overview <- data.frame(
  sheet = sheets,
  description = c(
    "Every protein and the module it was assigned to (grey = unassigned)",
    "Module sizes and how much eigengene variance subject identity still explains (ICC)",
    "GO:BP over-representation per module, universe = the measured proteome",
    "limma::fry: does each module move in each of the nine contrasts?",
    "Baseline (T1) eigengene vs each adaptation trait: in-sample r and LOSO q2",
    "ME ~ trait * timepoint + (1 | subject); the trait:timepoint term",
    "Module eigengene per sample, scored on raw abundance",
    "The soft-threshold sweep behind the chosen power"
  ),
  stringsAsFactors = FALSE
)

wb <- createWorkbook()
addWorksheet(wb, "overview")
writeData(wb, "overview", overview)
for (s in sheets) {
  addWorksheet(wb, s)
  writeData(wb, s, F04_AUDIT[[s]])
}
addWorksheet(wb, "metadata")
writeData(wb, "metadata", metadata)
saveWorkbook(wb, file.path(DAT_DIR, "WGCNA_source_data.xlsx"), overwrite = TRUE)

cat("F04 rebuilt: atlas, pathways, fry, phenotype, supplements, workbook\n")

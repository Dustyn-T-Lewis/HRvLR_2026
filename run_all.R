# Run stages in order, each in a separate R session (.R by Rscript, .qmd by quarto render);
# stop on the first failure.
# Rscript run_all.R          all stages
# Rscript run_all.R 01 04    stages starting with these prefixes; prerequisites must exist
# Logs are local to Scratch/runlogs/.

stages <- c(
  "01_Processing/01_Protein_Filtering",
  "01_Processing/02_Normalisation",
  "01_Processing/03_Quantification",
  "01_Processing/04_Gene_Sets",
  "01_Processing/05_Set_Scores",
  "01_Processing/06_Module_Building",
  "01_Processing/07_Module_Characterisation",
  "02_Groups/HRvLR",
  "02_Groups/Muscle_size",
  "02_Groups/Muscle_strength",
  "02_Groups/Proteome_T1",
  "02_Groups/Proteome_T2",
  "02_Groups/Proteome_T3",
  "02_Groups/Proteome_training",
  "02_Groups/Proteome_acute",
  "03_Group_comparisons/HRvLR",
  "03_Group_comparisons/Muscle_size",
  "03_Group_comparisons/Muscle_strength",
  "03_Group_comparisons/Proteome_T1",
  "03_Group_comparisons/Proteome_T2",
  "03_Group_comparisons/Proteome_T3",
  "03_Group_comparisons/Proteome_training",
  "03_Group_comparisons/Proteome_acute",
  # Reads all five proteome groupings, so it runs after them.
  "03_Group_comparisons/Proteome_clusters",
  "04_Continuous_no_group/01_Trait_Associations",
  "04_Continuous_no_group/02_Trait_Pathways",
  "04_Continuous_no_group/03_Profile_Clusters",
  # Figures read saved results only, so they run last.
  "05_Figures/Processing/S01_QC",
  "05_Figures/HRvLR/F1_Overview",
  "05_Figures/HRvLR/F2_Programmes",
  "05_Figures/HRvLR/F3_People",
  "05_Figures/HRvLR/S01_Volcano",
  "05_Figures/HRvLR/S02_NES_Scatter",
  "05_Figures/HRvLR/S03_Signatures",
  "05_Figures/Muscle_size/F1_Overview",
  "05_Figures/Muscle_strength/F1_Overview",
  "05_Figures/Proteome_T1/F1_Overview",
  "05_Figures/Proteome_T2/F1_Overview",
  "05_Figures/Proteome_T3/F1_Overview",
  "05_Figures/Proteome_training/F1_Overview",
  "05_Figures/Proteome_acute/F1_Overview",
  "05_Figures/Continuous/F1_Overview"
)
only <- commandArgs(trailingOnly = TRUE)
if (length(only)) stages <- stages[rowSums(outer(stages, only, startsWith)) > 0]
dir.create("Scratch/runlogs", recursive = TRUE, showWarnings = FALSE)
# A stage is a folder whose a_script/ holds one script: a .qmd is rendered into the folder's
# b_reports/, a .R is run with Rscript.
render_args <- function(stage) {
  script <- list.files(file.path(stage, "a_script"), pattern = "\\.(R|qmd)$", full.names = TRUE)
  stopifnot(length(script) == 1)
  if (endsWith(script, ".qmd")) c("render", script, "--output-dir", "../b_reports") else script
}
for (stage in stages) {
  args <- render_args(stage)
  log <- file.path("Scratch/runlogs", paste0(gsub("/", "_", stage, fixed = TRUE), ".log"))
  started <- Sys.time()
  # Called from the project root, quarto's R session loads the renv library before any package.
  program <- if (args[1] == "render") "quarto" else "Rscript"
  status <- system2(program, args, stdout = log, stderr = log)
  message(sprintf("%-32s %6.1f min", stage, as.numeric(Sys.time() - started, units = "mins")))
  if (status != 0) stop(stage, " failed; see ", log)
}

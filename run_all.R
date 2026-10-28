# Run stages in order, each in a separate R session (.R by Rscript, .qmd by quarto render);
# stop on the first failure.
# Rscript run_all.R          all stages
# Rscript run_all.R 01 04    stages starting with these prefixes; prerequisites must exist
# Logs are local to Scratch/runlogs/.

stages <- c(
  "01_Processing/01_Precursors",
  "01_Processing/02_Quantification",
  "01_Processing/03_Blood_Correction",
  "01_Processing/04_Gene_Sets",
  "01_Processing/05_Set_Scores",
  "01_Processing/06_Module_Building",
  "01_Processing/07_Module_Characterisation",
  "02_HRvLR/01_Protein_Contrasts",
  "02_HRvLR/02_Set_Contrasts",
  "02_HRvLR/03_Module_Contrasts",
  "02_HRvLR/04_Variance",
  "02_HRvLR/05_Module_Preservation",
  "02_HRvLR/06_Classification",
  "02_HRvLR/07_Group_Trait_Slopes",
  "02_HRvLR/08_Module_Detail",
  "03_Continuous/01_Trait_Associations",
  "03_Continuous/02_Trait_Pathways",
  "03_Continuous/03_Profile_Clusters",
  paste0("04_Figures/", c("F1", "S1", "F2", "S2", "F3", "S3", "F4", "S4", "F5", "S5", "Modules"))
)
only <- commandArgs(trailingOnly = TRUE)
if (length(only)) stages <- stages[rowSums(outer(stages, only, startsWith)) > 0]
dir.create("Scratch/runlogs", recursive = TRUE, showWarnings = FALSE)
# A stage is a folder whose a_script/ holds one script: a .qmd is rendered into the folder's
# b_reports/, a .R is run with Rscript. A figure's panel scripts sit one level down.
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

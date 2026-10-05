# Runs every substage in order, each in its own R session, and stops at the first failure.
#   Rscript run_all.R               all stages
#   Rscript run_all.R 03 05         only the stages whose folder starts with 03 or 05
# Logs go to .runlogs/<substage>.log.

stages <- c(
  "01_Preprocess/01_Filter", "01_Preprocess/02_Normalize", "01_Preprocess/03_Quantify",
  "02_Proteins/01_Design", "02_Proteins/02_Contrasts", "02_Proteins/03_Classify",
  "02_Proteins/04_Associate",
  "03_Pathways/01_Gene_Sets", "03_Pathways/02_Scores", "03_Pathways/03_Contrasts",
  "03_Pathways/04_Classify", "03_Pathways/05_Associate",
  "04_Modules/01_Build", "04_Modules/02_Contrasts", "04_Modules/03_Classify",
  "04_Modules/04_Associate", "04_Modules/05_Characterise", "04_Modules/06_Preserve",
  "05_Integration/01_Chance_Line", "05_Integration/02_Candidates", "05_Integration/03_Panel",
  "05_Integration/04_Continuous", "05_Integration/05_Clustering",
  "06_Figures/F01_Overview", "06_Figures/F02_Groups", "06_Figures/F03_Continuum",
  "06_Figures/S01_Volcano", "06_Figures/S02_NES_Scatter",
  "06_Figures/S03_Signatures", "06_Figures/S04_QC"
)
only <- commandArgs(trailingOnly = TRUE)
if (length(only)) stages <- stages[substr(stages, 1, 2) %in% only]
dir.create(".runlogs", showWarnings = FALSE)
for (stage in stages) {
  script <- list.files(file.path(stage, "a_script"), pattern = "\\.R$", full.names = TRUE)
  stopifnot(length(script) == 1)
  log <- file.path(".runlogs", paste0(basename(stage), ".log"))
  started <- Sys.time()
  status <- system2("Rscript", script, stdout = log, stderr = log)
  message(sprintf("%-32s %6.1f min", stage, as.numeric(Sys.time() - started, units = "mins")))
  if (status != 0) stop(stage, " failed; see ", log)
}

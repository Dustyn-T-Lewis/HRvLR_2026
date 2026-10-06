# Runs every substage in order, each in its own R session, and stops at the first failure.
#   Rscript run_all.R               all stages
#   Rscript run_all.R 01 04         only the arcs whose folder starts with 01 or 04
# Logs go to .runlogs/<arc>_<substage>.log.

stages <- c(
  "01_Trunk/01_Filter", "01_Trunk/02_Normalize", "01_Trunk/03_Quantify",
  "01_Trunk/04_Gene_Sets", "01_Trunk/05_Scores", "01_Trunk/06_Modules", "01_Trunk/07_Characterise",
  "04_Continuous_no_group/01_Associate/01_Proteins",
  "04_Continuous_no_group/01_Associate/02_Sets",
  "04_Continuous_no_group/01_Associate/03_Modules",
  "04_Continuous_no_group/02_Pathways", "04_Continuous_no_group/04_Profile_Clusters",
  "02_HRvLR_groups/01_Design", "02_HRvLR_groups/02_DEP", "02_HRvLR_groups/03_Pathways",
  "02_HRvLR_groups/04_Modules", "02_HRvLR_groups/05_Classify/01_Proteins",
  "02_HRvLR_groups/05_Classify/02_Sets", "02_HRvLR_groups/05_Classify/03_Modules",
  "02_HRvLR_groups/05_Classify/04_Panel", "02_HRvLR_groups/07_Preserve"
)
only <- commandArgs(trailingOnly = TRUE)
if (length(only)) stages <- stages[substr(stages, 1, 2) %in% only]
dir.create(".runlogs", showWarnings = FALSE)
for (stage in stages) {
  script <- list.files(file.path(stage, "a_script"), pattern = "\\.R$", full.names = TRUE)
  stopifnot(length(script) == 1)
  log <- file.path(".runlogs", paste0(gsub("/", "_", stage), ".log"))
  started <- Sys.time()
  status <- system2("Rscript", script, stdout = log, stderr = log)
  message(sprintf("%-32s %6.1f min", stage, as.numeric(Sys.time() - started, units = "mins")))
  if (status != 0) stop(stage, " failed; see ", log)
}

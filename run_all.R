# Runs every step in order, each in its own R session, and stops at the first failure. A full run
# then rewrites the step READMEs and pipeline.png.
#   Rscript run_all.R               all stages
#   Rscript run_all.R 01 04         only the folders starting 01 or 04
# Logs go to .runlogs/, one per step.

stages <- c(
  "01_Processing/01_Protein_Filtering",
  "01_Processing/02_Normalisation",
  "01_Processing/03_Quantification",
  "01_Processing/04_Gene_Sets",
  "01_Processing/05_Set_Scores",
  "01_Processing/06_Module_Building",
  "01_Processing/07_Module_Characterisation",
  "02_HRvLR_groups/01_Design_Contrasts",
  "02_HRvLR_groups/02_Protein_Contrasts",
  "02_HRvLR_groups/03_Set_Contrasts",
  "02_HRvLR_groups/04_Module_Contrasts",
  "02_HRvLR_groups/05_Module_Preservation",
  "02_HRvLR_groups/06_Classification/01_Proteins",
  "02_HRvLR_groups/06_Classification/02_Sets",
  "02_HRvLR_groups/06_Classification/03_Modules",
  "02_HRvLR_groups/07_Trait_Associations",
  "03_Cluster_groups/Muscle_size/01_Clustering",
  "03_Cluster_groups/Muscle_size/02_Design_Contrasts",
  "03_Cluster_groups/Muscle_size/03_Protein_Contrasts",
  "03_Cluster_groups/Muscle_size/04_Set_Contrasts",
  "03_Cluster_groups/Muscle_size/05_Module_Contrasts",
  "03_Cluster_groups/Muscle_size/06_Classification/01_Proteins",
  "03_Cluster_groups/Muscle_size/06_Classification/02_Sets",
  "03_Cluster_groups/Muscle_size/06_Classification/03_Modules",
  "03_Cluster_groups/Muscle_size/07_Trait_Associations",
  "03_Cluster_groups/Muscle_size/Figures/F1_Overview",
  "03_Cluster_groups/Muscle_strength/01_Clustering",
  "03_Cluster_groups/Muscle_strength/02_Design_Contrasts",
  "03_Cluster_groups/Muscle_strength/03_Protein_Contrasts",
  "03_Cluster_groups/Muscle_strength/04_Set_Contrasts",
  "03_Cluster_groups/Muscle_strength/05_Module_Contrasts",
  "03_Cluster_groups/Muscle_strength/06_Classification/01_Proteins",
  "03_Cluster_groups/Muscle_strength/06_Classification/02_Sets",
  "03_Cluster_groups/Muscle_strength/06_Classification/03_Modules",
  "03_Cluster_groups/Muscle_strength/07_Trait_Associations",
  "03_Cluster_groups/Muscle_strength/Figures/F1_Overview",
  "03_Cluster_groups/Proteome/01_Clustering",
  "03_Cluster_groups/Proteome/02_Design_Contrasts",
  "03_Cluster_groups/Proteome/03_Protein_Contrasts",
  "03_Cluster_groups/Proteome/04_Set_Contrasts",
  "03_Cluster_groups/Proteome/05_Module_Contrasts",
  "03_Cluster_groups/Proteome/06_Classification/01_Proteins",
  "03_Cluster_groups/Proteome/06_Classification/02_Sets",
  "03_Cluster_groups/Proteome/06_Classification/03_Modules",
  "03_Cluster_groups/Proteome/07_Trait_Associations",
  "03_Cluster_groups/Proteome/08_Cluster_Phenotype",
  "03_Cluster_groups/Proteome/Figures/F1_Overview",
  "04_Continuous_no_group/01_Trait_Associations/01_Proteins",
  "04_Continuous_no_group/01_Trait_Associations/02_Sets",
  "04_Continuous_no_group/01_Trait_Associations/03_Modules",
  "04_Continuous_no_group/02_Trait_Pathways",
  "04_Continuous_no_group/03_Profile_Clusters",
  "04_Continuous_no_group/Figures/F1_Overview",
  # Figures last: F3 reads 04_Continuous_no_group.
  "02_HRvLR_groups/Figures/F1_Overview",
  "02_HRvLR_groups/Figures/F2_Programmes",
  "02_HRvLR_groups/Figures/F3_People",
  "02_HRvLR_groups/Figures/S01_Volcano",
  "02_HRvLR_groups/Figures/S02_NES_Scatter",
  "02_HRvLR_groups/Figures/S03_Signatures",
  "01_Processing/Figures/S01_QC"
)
only <- commandArgs(trailingOnly = TRUE)
if (length(only)) stages <- stages[substr(stages, 1, 2) %in% only]
dir.create(".runlogs", showWarnings = FALSE)
run <- function(script, name) {
  log <- file.path(".runlogs", paste0(gsub("/", "_", name), ".log"))
  started <- Sys.time()
  status <- system2("Rscript", script, stdout = log, stderr = log)
  message(sprintf("%-32s %6.1f min", name, as.numeric(Sys.time() - started, units = "mins")))
  if (status != 0) stop(name, " failed; see ", log)
}
for (stage in stages) {
  script <- list.files(file.path(stage, "a_script"), pattern = "\\.R$", full.names = TRUE)
  stopifnot(length(script) == 1)
  run(script, stage)
}
# After a full run, the step READMEs and pipeline.png are redrawn from the new outputs.
if (!length(only)) {
  run("docs/readmes.R", "docs/readmes")
  run("docs/pipeline.R", "docs/pipeline")
}

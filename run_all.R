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
  "02_Groups/a_script/study_groups.qmd@HRvLR",
  "02_Groups/a_script/clustering.qmd@Muscle_size",
  "02_Groups/a_script/clustering.qmd@Muscle_strength",
  "02_Groups/a_script/clustering.qmd@Proteome",
  "03_Group_comparisons/a_script/01_design_contrasts.qmd@HRvLR",
  "03_Group_comparisons/a_script/02_protein_contrasts.qmd@HRvLR",
  "03_Group_comparisons/HRvLR/03_Set_Contrasts",
  "03_Group_comparisons/HRvLR/04_Module_Contrasts",
  "03_Group_comparisons/HRvLR/05_Classification/01_Proteins",
  "03_Group_comparisons/HRvLR/05_Classification/02_Sets",
  "03_Group_comparisons/HRvLR/05_Classification/03_Modules",
  "03_Group_comparisons/HRvLR/06_Trait_Associations",
  "03_Group_comparisons/HRvLR/07_Module_Preservation",
  "03_Group_comparisons/a_script/01_design_contrasts.qmd@Muscle_size",
  "03_Group_comparisons/a_script/02_protein_contrasts.qmd@Muscle_size",
  "03_Group_comparisons/Muscle_size/03_Set_Contrasts",
  "03_Group_comparisons/Muscle_size/04_Module_Contrasts",
  "03_Group_comparisons/Muscle_size/05_Classification/01_Proteins",
  "03_Group_comparisons/Muscle_size/05_Classification/02_Sets",
  "03_Group_comparisons/Muscle_size/05_Classification/03_Modules",
  "03_Group_comparisons/Muscle_size/06_Trait_Associations",
  "03_Group_comparisons/a_script/01_design_contrasts.qmd@Muscle_strength",
  "03_Group_comparisons/a_script/02_protein_contrasts.qmd@Muscle_strength",
  "03_Group_comparisons/Muscle_strength/03_Set_Contrasts",
  "03_Group_comparisons/Muscle_strength/04_Module_Contrasts",
  "03_Group_comparisons/Muscle_strength/05_Classification/01_Proteins",
  "03_Group_comparisons/Muscle_strength/05_Classification/02_Sets",
  "03_Group_comparisons/Muscle_strength/05_Classification/03_Modules",
  "03_Group_comparisons/Muscle_strength/06_Trait_Associations",
  "03_Group_comparisons/a_script/01_design_contrasts.qmd@Proteome/T1",
  "03_Group_comparisons/a_script/01_design_contrasts.qmd@Proteome/T2",
  "03_Group_comparisons/a_script/01_design_contrasts.qmd@Proteome/T3",
  "03_Group_comparisons/a_script/01_design_contrasts.qmd@Proteome/training",
  "03_Group_comparisons/a_script/01_design_contrasts.qmd@Proteome/acute",
  "03_Group_comparisons/a_script/02_protein_contrasts.qmd@Proteome/T1",
  "03_Group_comparisons/a_script/02_protein_contrasts.qmd@Proteome/T2",
  "03_Group_comparisons/a_script/02_protein_contrasts.qmd@Proteome/T3",
  "03_Group_comparisons/a_script/02_protein_contrasts.qmd@Proteome/training",
  "03_Group_comparisons/a_script/02_protein_contrasts.qmd@Proteome/acute",
  "03_Group_comparisons/Proteome/03_Set_Contrasts",
  "03_Group_comparisons/Proteome/04_Module_Contrasts",
  "03_Group_comparisons/Proteome/05_Classification/01_Proteins",
  "03_Group_comparisons/Proteome/05_Classification/02_Sets",
  "03_Group_comparisons/Proteome/05_Classification/03_Modules",
  "03_Group_comparisons/Proteome/06_Trait_Associations",
  "03_Group_comparisons/Proteome/07_Cluster_Phenotype",
  "04_Continuous_no_group/01_Trait_Associations/01_Proteins",
  "04_Continuous_no_group/01_Trait_Associations/02_Sets",
  "04_Continuous_no_group/01_Trait_Associations/03_Modules",
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
  "05_Figures/Proteome/F1_Overview",
  "05_Figures/Continuous/F1_Overview"
)
only <- commandArgs(trailingOnly = TRUE)
if (length(only)) stages <- stages[rowSums(outer(stages, only, startsWith)) > 0]
dir.create("Scratch/runlogs", recursive = TRUE, showWarnings = FALSE)
# A stage is a step folder holding one script in a_script/, or "<walkthrough>.qmd@<grouping>",
# optionally "@<grouping>/<input>": one walkthrough rendered with params$grouping (and $input)
# set. A numbered walkthrough (02_protein_contrasts.qmd) writes its report into its step folder
# (02_Protein_Contrasts) under the grouping; an unnumbered one into the grouping folder itself.
render_args <- function(stage) {
  if (grepl("@", stage)) {
    script <- sub("@.*", "", stage)
    target <- strsplit(sub(".*@", "", stage), "/")[[1]]
    stem <- tools::file_path_sans_ext(basename(script))
    step <- if (grepl("^[0-9]", stem)) gsub("_([a-z])", "_\\U\\1", stem, perl = TRUE)
    params <- c("-P", paste0("grouping:", target[1]))
    if (length(target) > 1) params <- c(params, "-P", paste0("input:", target[2]))
    out <- do.call(file.path, as.list(c("..", target[1], step, target[-1], "b_reports")))
    return(c("render", script, params, "--output-dir", out))
  }
  script <- list.files(file.path(stage, "a_script"), pattern = "\\.(R|qmd)$", full.names = TRUE)
  stopifnot(length(script) == 1)
  if (endsWith(script, ".qmd")) c("render", script, "--output-dir", "../b_reports") else script
}
for (stage in stages) {
  args <- render_args(stage)
  log <- file.path("Scratch/runlogs", paste0(gsub("[/@]", "_", stage), ".log"))
  started <- Sys.time()
  # Called from the project root, quarto's R session loads the renv library before any package.
  program <- if (args[1] == "render") "quarto" else "Rscript"
  status <- system2(program, args, stdout = log, stderr = log)
  message(sprintf("%-32s %6.1f min", stage, as.numeric(Sys.time() - started, units = "mins")))
  if (status != 0) stop(stage, " failed; see ", log)
}

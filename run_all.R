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
  "03_Group_comparisons/a_script/01_design_contrasts.qmd@HRvLR",
  "03_Group_comparisons/a_script/02_protein_contrasts.qmd@HRvLR",
  "03_Group_comparisons/a_script/03_set_contrasts.qmd@HRvLR",
  "03_Group_comparisons/a_script/04_module_contrasts.qmd@HRvLR",
  "03_Group_comparisons/a_script/05_classification.qmd@HRvLR",
  "03_Group_comparisons/a_script/06_trait_associations.qmd@HRvLR",
  "03_Group_comparisons/a_script/07_module_preservation.qmd@HRvLR",
  "03_Group_comparisons/a_script/01_design_contrasts.qmd@Muscle_size",
  "03_Group_comparisons/a_script/02_protein_contrasts.qmd@Muscle_size",
  "03_Group_comparisons/a_script/03_set_contrasts.qmd@Muscle_size",
  "03_Group_comparisons/a_script/04_module_contrasts.qmd@Muscle_size",
  "03_Group_comparisons/a_script/05_classification.qmd@Muscle_size",
  "03_Group_comparisons/a_script/06_trait_associations.qmd@Muscle_size",
  "03_Group_comparisons/a_script/01_design_contrasts.qmd@Muscle_strength",
  "03_Group_comparisons/a_script/02_protein_contrasts.qmd@Muscle_strength",
  "03_Group_comparisons/a_script/03_set_contrasts.qmd@Muscle_strength",
  "03_Group_comparisons/a_script/04_module_contrasts.qmd@Muscle_strength",
  "03_Group_comparisons/a_script/05_classification.qmd@Muscle_strength",
  "03_Group_comparisons/a_script/06_trait_associations.qmd@Muscle_strength",
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
  "03_Group_comparisons/a_script/03_set_contrasts.qmd@Proteome/T1",
  "03_Group_comparisons/a_script/03_set_contrasts.qmd@Proteome/T2",
  "03_Group_comparisons/a_script/03_set_contrasts.qmd@Proteome/T3",
  "03_Group_comparisons/a_script/03_set_contrasts.qmd@Proteome/training",
  "03_Group_comparisons/a_script/03_set_contrasts.qmd@Proteome/acute",
  "03_Group_comparisons/a_script/04_module_contrasts.qmd@Proteome/T1",
  "03_Group_comparisons/a_script/04_module_contrasts.qmd@Proteome/T2",
  "03_Group_comparisons/a_script/04_module_contrasts.qmd@Proteome/T3",
  "03_Group_comparisons/a_script/04_module_contrasts.qmd@Proteome/training",
  "03_Group_comparisons/a_script/04_module_contrasts.qmd@Proteome/acute",
  "03_Group_comparisons/a_script/05_classification.qmd@Proteome/T1",
  "03_Group_comparisons/a_script/05_classification.qmd@Proteome/T2",
  "03_Group_comparisons/a_script/05_classification.qmd@Proteome/T3",
  "03_Group_comparisons/a_script/05_classification.qmd@Proteome/training",
  "03_Group_comparisons/a_script/05_classification.qmd@Proteome/acute",
  "03_Group_comparisons/a_script/06_trait_associations.qmd@Proteome/T1",
  "03_Group_comparisons/a_script/06_trait_associations.qmd@Proteome/T2",
  "03_Group_comparisons/a_script/06_trait_associations.qmd@Proteome/T3",
  "03_Group_comparisons/a_script/06_trait_associations.qmd@Proteome/training",
  "03_Group_comparisons/a_script/06_trait_associations.qmd@Proteome/acute",
  "03_Group_comparisons/a_script/07_cluster_phenotype.qmd@Proteome",
  "04_Continuous_no_group/01_Trait_Associations",
  "04_Continuous_no_group/02_Trait_Pathways",
  "04_Continuous_no_group/03_Profile_Clusters",
  # Figures read saved results only, so they run last.
  "05_Figures/Processing/S01_QC",
  "05_Figures/a_script/F1_Overview.qmd@HRvLR",
  "05_Figures/HRvLR/F2_Programmes",
  "05_Figures/HRvLR/F3_People",
  "05_Figures/HRvLR/S01_Volcano",
  "05_Figures/HRvLR/S02_NES_Scatter",
  "05_Figures/HRvLR/S03_Signatures",
  "05_Figures/a_script/F1_Overview.qmd@Muscle_size",
  "05_Figures/a_script/F1_Overview.qmd@Muscle_strength",
  "05_Figures/a_script/F1_Overview.qmd@Proteome",
  "05_Figures/Continuous/F1_Overview"
)
only <- commandArgs(trailingOnly = TRUE)
if (length(only)) stages <- stages[rowSums(outer(stages, only, startsWith)) > 0]
dir.create("Scratch/runlogs", recursive = TRUE, showWarnings = FALSE)
# A stage is a step folder holding one script in a_script/, or "<walkthrough>.qmd@<grouping>",
# optionally "@<grouping>/<input>": one walkthrough rendered with params$grouping (and $input)
# set. A numbered walkthrough (02_protein_contrasts.qmd) writes its report into its step folder
# (02_Protein_Contrasts) under the grouping, a capitalised one (F1_Overview.qmd) into the folder of
# its own name, and any other into the grouping folder itself.
render_args <- function(stage) {
  if (grepl("@", stage)) {
    script <- sub("@.*", "", stage)
    target <- strsplit(sub(".*@", "", stage), "/")[[1]]
    stem <- tools::file_path_sans_ext(basename(script))
    step <- if (grepl("^[0-9]", stem)) {
      gsub("_([a-z])", "_\\U\\1", stem, perl = TRUE)
    } else if (grepl("^[A-Z]", stem)) {
      stem
    }
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

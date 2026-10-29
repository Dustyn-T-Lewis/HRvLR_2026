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
  "03_Continuous/04_Elastic_Net",
  paste0("04_Figures/", c("F1", "F2", "F3", "F4", "F5"))
)
only <- commandArgs(trailingOnly = TRUE)
if (length(only)) stages <- stages[rowSums(outer(stages, only, startsWith)) > 0]
dir.create("Scratch/runlogs", recursive = TRUE, showWarnings = FALSE)
# A step's a_script/ holds one .qmd, rendered into its b_reports/. A figure's a_script/ holds its
# panel scripts and <F>_composite.R, then supp/ with each supplement's composite; every composite
# runs with Rscript, the figure's first.
scripts_of <- function(stage) {
  qmd <- list.files(file.path(stage, "a_script"), "\\.qmd$", full.names = TRUE)
  if (length(qmd)) {
    return(qmd)
  }
  c(
    list.files(file.path(stage, "a_script"), "_composite\\.R$", full.names = TRUE),
    list.files(file.path(stage, "a_script", "supp"), "_composite\\.R$", full.names = TRUE)
  )
}
for (stage in stages) {
  for (script in scripts_of(stage)) {
    name <- sub("\\.(R|qmd)$", "", basename(script))
    log <- file.path(
      "Scratch/runlogs", paste0(gsub("/", "_", stage, fixed = TRUE), "_", name, ".log")
    )
    started <- Sys.time()
    # Called from the project root, quarto's R session loads the renv library before any package.
    status <- if (endsWith(script, ".qmd")) {
      system2("quarto", c("render", script, "--output-dir", "../b_reports"),
        stdout = log, stderr = log
      )
    } else {
      system2("Rscript", script, stdout = log, stderr = log)
    }
    minutes <- as.numeric(Sys.time() - started, units = "mins")
    message(sprintf("%-44s %6.1f min", paste(stage, name), minutes))
    if (status != 0) stop(stage, " failed; see ", log)
  }
}

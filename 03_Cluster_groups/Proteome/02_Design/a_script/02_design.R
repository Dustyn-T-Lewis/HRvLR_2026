# Run once per proteome clustering (T1, T2, T3, training change, acute change) on that
# clustering's groups C1 and C2: a copy of 02_HRvLR_groups/01_Design. Each input writes to its own
# folder; stack_inputs() then joins them into c_data/all_inputs.xlsx with a validity tag.
# Six cell means and the blood index with subject blocked, the nine contrasts and their roles, and
# the within-subject correlation limpa's fit estimates. Every later limma fit reads design.rds.

suppressPackageStartupMessages({
  library(here)
  library(limma)
  library(dplyr)
  library(tibble)
  library(writexl)
})
source(here("R", "helpers.R"))
source(here("03_Cluster_groups", "Proteome", "R", "inputs.R"))

run_input <- function(input) {
  inputs <- c(proteins = "03_Cluster_groups/Proteome/01_Cluster/c_data/inputs/{input}/proteins.rds")
  inputs <- sub("{input}", input, inputs, fixed = TRUE)
  paths <- vapply(inputs, here, character(1))
  stopifnot(file.exists(paths))
  quantified <- readRDS(paths[["proteins"]])
  proteins <- quantified$proteins
  meta <- quantified$metadata |>
    mutate(
      arm = factor(arm, levels = c("C1", "C2")),
      timepoint = factor(timepoint, levels = c("T1", "T2", "T3"))
    )
  stopifnot(
    identical(meta$sample_id, colnames(proteins$E)), !anyNA(meta$group),
    !any(is.na(meta$subject) | meta$subject == "")
  )

  # C1 and C2 are different people, so subject is blocked rather than fixed: a fixed subject term
  # would absorb the between-arm contrasts. Blood enters as a covariate.
  design <- model.matrix(~ 0 + group + blood, meta)
  colnames(design) <- sub("^group", "", colnames(design))
  rownames(design) <- meta$sample_id

  # Rank deficiency would make the contrasts inestimable and surface as silent NAs.
  stopifnot(
    qr(design)$rank == ncol(design),
    identical(make.names(colnames(design)), colnames(design))
  )

  contrast_vector <- c(
    "Training_C1 = C1_T2 - C1_T1",
    "Training_C2 = C2_T2 - C2_T1",
    "Acute_C1 = C1_T3 - C1_T2",
    "Acute_C2 = C2_T3 - C2_T2",
    "Baseline_C1vC2 = C1_T1 - C2_T1",
    "Trained_C1vC2 = C1_T2 - C2_T2",
    "Acute_C1vC2 = C1_T3 - C2_T3",
    "Training_Interaction = (C1_T2 - C1_T1) - (C2_T2 - C2_T1)",
    "Acute_Interaction = (C1_T3 - C1_T2) - (C2_T3 - C2_T2)"
  )
  contrast_matrix <- makeContrasts(contrasts = sub(".*= ", "", contrast_vector), levels = design)
  contrast_names <- sub(" =.*", "", contrast_vector)
  colnames(contrast_matrix) <- contrast_names

  # Baseline_C1vC2 is the floor, not a negative control: the arm label comes from the outcome.
  roles <- tibble(
    contrast = contrast_names,
    role = case_when(
      contrast == "Training_Interaction" ~ "primary",
      contrast == "Acute_Interaction" ~ "secondary",
      contrast == "Baseline_C1vC2" ~ "floor",
      TRUE ~ "descriptive"
    )
  )

  strata <- tibble(
    level = "within subject, from limpa's dpcDE fit with subject blocked",
    correlation = limpa::dpcDE(proteins, design, block = meta$subject, plot = FALSE)$correlation
  )
  # Downstream steps read the object proteoDA used to write: data, annotation, metadata, design.
  dal <- list(
    data = proteins$E, annotation = quantified$annotation, metadata = meta,
    design = list(design_matrix = design, contrast_matrix = contrast_matrix)
  )

  out <- here("03_Cluster_groups", "Proteome", "02_Design", input, "c_data")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  saveRDS(
    list(dal = dal, contrasts = contrast_vector, roles = roles, floor = "Baseline_C1vC2"),
    file.path(out, "design.rds"),
    compress = "xz"
  )
  sheets <- list(
    contrasts = tibble(definition = contrast_vector, contrast = contrast_names) |>
      left_join(roles, by = "contrast"),
    correlation_strata = strata
  )
  holds <- c(
    "The nine contrasts, their definitions and roles.",
    "Within-subject correlation from limpa's fit, blocked by subject."
  )
  write_workbook(sheets, holds, file.path(out, "01_design.xlsx"), inputs, paths)
}
purrr::walk(inputs_used, run_input)
stack_inputs(here("03_Cluster_groups", "Proteome", "02_Design"))

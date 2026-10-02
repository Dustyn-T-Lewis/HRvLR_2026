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

inputs <- c(proteins = "01_Preprocess/03_Quantification/c_data/proteins.rds")
paths <- vapply(inputs, here, character(1))
stopifnot(file.exists(paths))
quantified <- readRDS(paths[["proteins"]])
proteins <- quantified$proteins
meta <- quantified$metadata |>
  mutate(
    arm = factor(arm, levels = c("HR", "LR")),
    timepoint = factor(timepoint, levels = c("T1", "T2", "T3"))
  )
stopifnot(
  identical(meta$sample_id, colnames(proteins$E)), !anyNA(meta$group),
  !any(is.na(meta$subject) | meta$subject == "")
)

# HR and LR are different people, so subject is blocked rather than fixed: a fixed subject term
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
  "Training_HR = HR_T2 - HR_T1",
  "Training_LR = LR_T2 - LR_T1",
  "Acute_HR = HR_T3 - HR_T2",
  "Acute_LR = LR_T3 - LR_T2",
  "Baseline_HRvLR = HR_T1 - LR_T1",
  "Trained_HRvLR = HR_T2 - LR_T2",
  "Acute_HRvLR = HR_T3 - LR_T3",
  "Training_Interaction = (HR_T2 - HR_T1) - (LR_T2 - LR_T1)",
  "Acute_Interaction = (HR_T3 - HR_T2) - (LR_T3 - LR_T2)"
)
contrast_matrix <- makeContrasts(contrasts = sub(".*= ", "", contrast_vector), levels = design)
contrast_names <- sub(" =.*", "", contrast_vector)
colnames(contrast_matrix) <- contrast_names

# Baseline_HRvLR is the floor, not a negative control: the arm label comes from the outcome.
roles <- tibble(
  contrast = contrast_names,
  role = case_when(
    contrast == "Training_Interaction" ~ "primary",
    contrast == "Acute_Interaction" ~ "secondary",
    contrast == "Baseline_HRvLR" ~ "floor",
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

out <- here("02_Differential_Expression", "01_Design", "c_data")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
saveRDS(
  list(dal = dal, contrasts = contrast_vector, roles = roles, floor = "Baseline_HRvLR"),
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

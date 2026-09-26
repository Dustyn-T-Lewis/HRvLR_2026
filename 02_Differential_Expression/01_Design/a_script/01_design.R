# Six cell means with subject random, the nine contrasts and their roles, and the within-subject
# correlation. Every later limma fit reads design.rds.

suppressPackageStartupMessages({
  library(here)
  library(proteoDA)
  library(limma)
  library(dplyr)
  library(tibble)
  library(writexl)
})

inputs <- c(normalized = "01_Preprocess/02_Normalization/c_data/DAList_normalized.rds")
paths <- vapply(inputs, here, character(1))
stopifnot(file.exists(paths))
dal <- readRDS(paths[["normalized"]])
stopifnot(identical(dal$metadata$sample_id, colnames(dal$data)))

group_levels <- c("HR_T1", "HR_T2", "HR_T3", "LR_T1", "LR_T2", "LR_T3")
meta <- dal$metadata |>
  mutate(
    arm = factor(arm, levels = c("HR", "LR")),
    timepoint = factor(timepoint, levels = c("T1", "T2", "T3")),
    group = factor(group, levels = group_levels)
  )
rownames(meta) <- meta$sample_id
stopifnot(!anyNA(meta$group), !any(is.na(meta$subject) | meta$subject == ""))
dal$metadata <- meta

# HR and LR are different people, so subject is a random effect: a fixed subject term would absorb
# the between-arm contrasts.
dal <- add_design(dal, "~ 0 + group + (1 | subject)")
design <- dal$design$design_matrix

# Rank deficiency would make the contrasts inestimable and surface as silent NAs.
# add_contrasts() parses the contrast strings, so every column name must survive make.names().
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
dal <- add_contrasts(dal, contrasts_vector = contrast_vector)
contrast_names <- colnames(dal$design$contrast_matrix)

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
  level = "within subject, subject random (unweighted matrix correlation)",
  correlation = duplicateCorrelation(dal$data, design, block = meta$subject)$consensus.correlation
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
  correlation_strata = strata,
  input_manifest = tibble(
    input = names(inputs), path = unname(inputs), md5 = unname(tools::md5sum(paths))
  ),
  package_versions = sessioninfo::package_info("loaded", dependencies = FALSE) |>
    as_tibble() |>
    select(package, version = loadedversion, source)
)
read_me <- tibble(sheet = names(sheets), holds = c(
  "The nine contrasts, their definitions and roles.",
  "Within-subject correlation from duplicateCorrelation, blocked by subject.",
  "Files read, with md5.",
  "Packages loaded at run time."
))
write_xlsx(c(list(read_me = read_me), sheets), file.path(out, "01_design.xlsx"))

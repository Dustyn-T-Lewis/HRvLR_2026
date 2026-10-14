#' ---
#' title: "Design and contrasts"
#' params:
#'   grouping: HRvLR
#' format:
#'   typst:
#'     toc: true
#'     margin:
#'       x: 1.6cm
#' execute:
#'   echo: true
#'   warning: false
#' ---
#'
#' Every comparison between two groups in this repository uses one model and one set of nine
#' contrasts. This step builds them for one grouping and saves them as `01_design.rds`, which every
#' later step in `03_Group_comparisons` reads.
#'
#' ## Packages
#'
#' limpa's fit estimates how alike one person's biopsies are; limma builds the contrasts.

#| label: packages
#| message: false
library(here)
library(limma)
library(dplyr)
library(purrr)
library(tibble)
source(here("00_Setup", "Shared_functions", "analysis_reporting_helpers.R"))
params$grouping

#' ## Inputs
#'
#' `proteins.rds` holds limpa's values for every protein in every run, with the sample sheet.
#' `groups.rds` from `02_Groups` says who is in each group and which group comes first.

#| label: inputs
data_dir <- file.path("03_Group_comparisons", params$grouping, "c_data")
inputs <- c(
  proteins = "01_Processing/03_Quantification/c_data/proteins.rds",
  groups = file.path("02_Groups", params$grouping, "c_data/groups.rds")
)
paths <- map_chr(inputs, here)
stopifnot(file.exists(paths))
quantified <- readRDS(paths[["proteins"]])
grouping <- readRDS(paths[["groups"]])
arms <- grouping$arms
table(grouping$group_of)

#' Each run takes its person's group label. People the grouping leaves out (for example a person a
#' clustering could not place) drop out with all their runs. `group` combines label and time point,
#' for example `HR_T1`.

#| label: relabel
meta <- relabel_samples(quantified$metadata, grouping$group_of) |>
  mutate(
    arm = factor(arm, levels = arms),
    timepoint = factor(timepoint, levels = c("T1", "T2", "T3"))
  )
proteins <- quantified$proteins[, meta$sample_id]
stopifnot(
  identical(meta$sample_id, colnames(proteins$E)), !anyNA(meta$group),
  !any(is.na(meta$subject) | meta$subject == "")
)
count(meta, group)

#' ## The model
#'
#' A cell-means model: one mean per group-by-time cell (six cells) plus the blood index as a
#' covariate, so every comparison is made at equal blood content. The groups are different people,
#' so person enters as a blocking factor rather than a fixed effect: a fixed person term would
#' absorb every between-group difference. Blocking tells the fit how alike one person's biopsies
#' are.

#| label: design
design <- model.matrix(~ 0 + group + blood, meta)
colnames(design) <- sub("^group", "", colnames(design))
rownames(design) <- meta$sample_id
stopifnot(
  qr(design)$rank == ncol(design),
  identical(make.names(colnames(design)), colnames(design))
)
head(design, 3)

#' The rank check stops the script if any contrast could not be estimated; otherwise that would
#' surface later as silent `NA`s.
#'
#' ## Nine contrasts
#'
#' A contrast is one comparison of cell means. With the first group as A and the second as B:
#'
#' | contrast | compares | window |
#' |---|---|---|
#' | `Training_A`, `Training_B` | T2 − T1 within each group | over training |
#' | `Acute_A`, `Acute_B` | T3 − T2 within each group | after the bout |
#' | `Baseline_AvB` | A − B at T1 | before training |
#' | `Trained_AvB` | A − B at T2 | after training |
#' | `Acute_AvB` | A − B at T3 | after the bout |
#' | `Training_Interaction` | A's training change − B's | over training |
#' | `Acute_Interaction` | A's acute change − B's | after the bout |
#'
#' The training interaction is the primary contrast and the acute interaction the secondary.
#' `Baseline_AvB` is the floor, not a negative control: the groups come from an outcome or from the
#' data, so a difference before training shows how far apart the split puts people to begin with.

#| label: contrasts
contrast_vector <- stringr::str_replace_all(
  c(
    "Training_{a} = {a}_T2 - {a}_T1",
    "Training_{b} = {b}_T2 - {b}_T1",
    "Acute_{a} = {a}_T3 - {a}_T2",
    "Acute_{b} = {b}_T3 - {b}_T2",
    "Baseline_{a}v{b} = {a}_T1 - {b}_T1",
    "Trained_{a}v{b} = {a}_T2 - {b}_T2",
    "Acute_{a}v{b} = {a}_T3 - {b}_T3",
    "Training_Interaction = ({a}_T2 - {a}_T1) - ({b}_T2 - {b}_T1)",
    "Acute_Interaction = ({a}_T3 - {a}_T2) - ({b}_T3 - {b}_T2)"
  ),
  c("\\{a\\}" = arms[1], "\\{b\\}" = arms[2])
)
contrast_matrix <- makeContrasts(contrasts = sub(".*= ", "", contrast_vector), levels = design)
contrast_names <- sub(" =.*", "", contrast_vector)
colnames(contrast_matrix) <- contrast_names
floor <- paste0("Baseline_", arms[1], "v", arms[2])
roles <- tibble(
  contrast = contrast_names,
  role = case_when(
    contrast == "Training_Interaction" ~ "primary",
    contrast == "Acute_Interaction" ~ "secondary",
    contrast == floor ~ "floor",
    TRUE ~ "descriptive"
  )
)
roles

#' ## How alike one person's biopsies are
#'
#' limpa's fit (`dpcDE()` with person as the block) estimates the consensus correlation between
#' biopsies from the same person across all proteins. Later protein-level fits reuse the same design
#' and block.

#| label: correlation
strata <- tibble(
  level = "within subject, from limpa's dpcDE fit with subject blocked",
  correlation = limpa::dpcDE(proteins, design, block = meta$subject, plot = FALSE)$correlation
)
strata

#' ## Validation

#| label: validate
stopifnot(
  ncol(contrast_matrix) == 9,
  all(colSums(abs(contrast_matrix[rownames(contrast_matrix) != "blood", ])) > 0),
  all(contrast_matrix["blood", ] == 0),
  setequal(unique(meta$arm), arms)
)

#' ## Summary
#'
#' `r params$grouping`: `r length(unique(meta$subject))` people,
#' `r nrow(meta)` runs, groups `r paste(arms, collapse = " and ")`. Within-person correlation
#' `r round(strata$correlation, 3)`.
#'
#' ## Save
#'
#' `01_design.rds` carries the protein matrix in the object later steps read (`dal`: data,
#' annotation, metadata, design), the contrasts and their roles, and the grouping.

#| label: save
dal <- list(
  data = proteins$E, annotation = quantified$annotation, metadata = meta,
  design = list(design_matrix = design, contrast_matrix = contrast_matrix)
)
out <- here(data_dir)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
saveRDS(
  list(
    dal = dal, contrasts = contrast_vector, roles = roles, floor = floor,
    group_of = grouping$group_of, arms = arms, defining_traits = grouping$defining_traits
  ),
  file.path(out, "01_design.rds"),
  compress = "xz"
)
write_workbook(
  list(
    contrasts = tibble(definition = contrast_vector, contrast = contrast_names) |>
      left_join(roles, by = "contrast"),
    correlation_strata = strata
  ),
  c(
    "The nine contrasts, their definitions and roles.",
    "Within-subject correlation from limpa's fit, blocked by subject."
  ),
  file.path(out, "01_design_contrasts.xlsx"), inputs, paths
)

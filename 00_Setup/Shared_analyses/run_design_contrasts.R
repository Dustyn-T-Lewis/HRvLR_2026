# Six cell means and the blood index with subject blocked, the nine contrasts and their roles, and
# the within-subject correlation limpa's fit estimates. Every later step of an arc reads design.rds:
# the design, the two groups and who is in each.
#
# groups is a clustering step's groups.rds (group_of, arms, defining_traits), or NULL for the
# study's own HR and LR. The first arm is the one each contrast puts first.
run_design_contrasts <- function(stage, groups = NULL) {
  inputs <- c(proteins = "01_Processing/03_Quantification/c_data/proteins.rds")
  if (!is.null(groups)) inputs <- c(inputs, groups = groups)
  paths <- map_chr(inputs, here)
  stopifnot(file.exists(paths))
  quantified <- readRDS(paths[["proteins"]])

  grouping <- if (is.null(groups)) {
    list(
      group_of = with(distinct(quantified$metadata, subject, arm), set_names(arm, subject)),
      arms = c("HR", "LR"), defining_traits = "comp_hypertrophy"
    )
  } else {
    readRDS(paths[["groups"]])
  }
  arms <- grouping$arms
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

  # The groups are different people, so subject is blocked rather than fixed: a fixed subject term
  # would absorb the between-arm contrasts. Blood enters as a covariate.
  design <- model.matrix(~ 0 + group + blood, meta)
  colnames(design) <- sub("^group", "", colnames(design))
  rownames(design) <- meta$sample_id

  # Rank deficiency would make the contrasts inestimable and surface as silent NAs.
  stopifnot(
    qr(design)$rank == ncol(design),
    identical(make.names(colnames(design)), colnames(design))
  )

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

  # The baseline contrast is the floor, not a negative control: the groups come from the outcome
  # or from the data.
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

  strata <- tibble(
    level = "within subject, from limpa's dpcDE fit with subject blocked",
    correlation = limpa::dpcDE(proteins, design, block = meta$subject, plot = FALSE)$correlation
  )
  # Later steps read the object proteoDA used to write: data, annotation, metadata, design.
  dal <- list(
    data = proteins$E, annotation = quantified$annotation, metadata = meta,
    design = list(design_matrix = design, contrast_matrix = contrast_matrix)
  )

  out <- here(stage, "c_data")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  saveRDS(
    list(
      dal = dal, contrasts = contrast_vector, roles = roles, floor = floor,
      group_of = grouping$group_of, arms = arms, defining_traits = grouping$defining_traits
    ),
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
  write_workbook(sheets, holds, file.path(out, "design_contrasts.xlsx"), inputs, paths)
}

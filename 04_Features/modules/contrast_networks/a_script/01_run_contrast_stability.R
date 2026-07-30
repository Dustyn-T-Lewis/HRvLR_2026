#!/usr/bin/env Rscript
# Is a contrast-specific co-expression network usable at this cohort size?
#
# The main network is built on 45 within-subject-centred samples and only 2 of
# its 11 modules survive leave-one-subject-out rebuilding. A training or acute
# network carries one value per subject instead, so it has roughly a third of
# the observations. This measures what that costs, using the same subsampling
# test, so the question is closed with a number rather than an argument.
#
# fit_modules() is deliberately not reused: it centres within subject, and a
# delta matrix has one column per subject, so centring would zero it entirely.

pacman::p_load(here, dplyr, openxlsx, WGCNA)
source(here("functions", "shared_wgcna.R"))
source(here("functions", "loso_wgcna_refit.R"))

set.seed(42)
disableWGCNAThreads()

imputed <- readRDS(here(
  "02_Normalization", "imputation", "c_data", "DAList_imputed_missforest.rds"
))
abund <- as.matrix(imputed$data)
rownames(abund) <- imputed$annotation$gene
meta <- imputed$metadata

# One column per subject: the change between two timepoints. A subject missing
# either timepoint cannot contribute.
contrast_matrix <- function(from, to) {
  subj <- unique(meta$Subject_ID)
  cols <- lapply(subj, function(s) {
    a <- meta$Col_ID[meta$Subject_ID == s & meta$Timepoint == from]
    b <- meta$Col_ID[meta$Subject_ID == s & meta$Timepoint == to]
    if (!length(a) || !length(b)) {
      return(NULL)
    }
    abund[, b] - abund[, a]
  })
  names(cols) <- subj
  cols <- Filter(Negate(is.null), cols)
  do.call(cbind, cols)
}

# Same parameters as the cohort network, minus the centring that a one-column
# -per-subject matrix cannot support.
fit_contrast_modules <- function(expr_t) {
  pw <- choose_power(expr_t)
  net <- do.call(WGCNA::blockwiseModules, c(
    list(
      expr_t,
      power = pw$power,
      networkType = "signed", TOMType = "signed",
      deepSplit = 2, minModuleSize = 30, mergeCutHeight = 0.15,
      pamStage = TRUE, pamRespectsDendro = TRUE,
      maxBlockSize = 2500, numericLabels = FALSE,
      randomSeed = 54321, verbose = 0
    ),
    BICOR_ARGS
  ))
  list(colors = setNames(net$colors, colnames(expr_t)), power = pw$power)
}

stability_of <- function(label, from, to) {
  x <- contrast_matrix(from, to)
  message(sprintf("%s: %d proteins x %d subjects", label, nrow(x), ncol(x)))
  full <- fit_contrast_modules(t(x))
  mods <- setdiff(unique(full$colors), "grey")
  message(sprintf(
    "  power %d, %d modules, %d of %d proteins in grey",
    full$power, length(mods), sum(full$colors == "grey"), length(full$colors)
  ))
  if (!length(mods)) {
    return(list(
      summary = data.frame(
        contrast = label, n_subjects = ncol(x), power = full$power,
        n_modules = 0L, n_stable = 0L, median_jaccard = NA_real_
      ),
      stability = NULL
    ))
  }
  folds <- lapply(colnames(x), function(s) {
    keep <- colnames(x) != s
    net <- fit_contrast_modules(t(x[, keep, drop = FALSE]))
    match_modules(net$colors, full$colors)
  })
  j <- bind_rows(folds) |>
    group_by(.data$full) |>
    summarise(
      mean_jaccard = mean(.data$jaccard), min_jaccard = min(.data$jaccard),
      n_folds = dplyr::n(), n_failed = sum(.data$jaccard < MIN_JACCARD),
      .groups = "drop"
    ) |>
    rename(module = "full") |>
    mutate(contrast = label, .before = 1)
  list(
    summary = data.frame(
      contrast = label, n_subjects = ncol(x), power = full$power,
      n_modules = length(mods), n_stable = sum(j$n_failed == 0),
      median_jaccard = median(j$mean_jaccard)
    ),
    stability = j
  )
}

results <- list(
  training = stability_of("training (T2-T1)", "T1", "T2"),
  acute = stability_of("acute (T3-T2)", "T2", "T3")
)

summary_tbl <- bind_rows(lapply(results, `[[`, "summary"))
stability_tbl <- bind_rows(lapply(results, `[[`, "stability"))

out <- here("04_Features", "modules", "contrast_networks", "c_data")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
wb <- createWorkbook()
addWorksheet(wb, "summary")
writeData(wb, "summary", summary_tbl)
if (nrow(stability_tbl)) {
  addWorksheet(wb, "module_stability")
  writeData(wb, "module_stability", stability_tbl)
}
saveWorkbook(wb, file.path(out, "contrast_stability.xlsx"), overwrite = TRUE)

print(summary_tbl, row.names = FALSE)

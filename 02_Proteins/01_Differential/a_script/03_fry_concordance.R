# Does the HR response signature also move in LR?
#
# The nine contrasts ask whether the arms differ. This asks the complement:
# take the proteins that responded in one arm, and test whether they move as a
# set in the other arm's ranking. A significant result is concordance, not
# difference, and concordance is what would explain an empty interaction.
#
# fry rather than fgsea. It is a rotation test, so it takes the design, the
# subject block and the consensus correlation, and it does not assume proteins
# are exchangeable. Two of this project's own findings depend on that
# distinction: fgsea's gene-permutation null scored random labels as well as
# real ones on a between-group contrast, and V1 replaced an fgsea module test
# with fry after the fgsea version returned padj = 7e-33 on a circular design.
#
# Only two pairings are testable without a calibration argument. A set and a
# ranking must share no cell mean, or the set partly defines its own ranking:
#
#   Training_HR (HR_T2 - HR_T1) vs Training_LR (LR_T2 - LR_T1)
#   Acute_HR    (HR_T3 - HR_T2) vs Acute_LR    (LR_T3 - LR_T2)
#
# Baseline_HRvLR shares HR_T1 with Training_HR and is excluded for that reason.

pacman::p_load(here, dplyr, tibble, readr, limma, openxlsx)

source(here("functions", "contrasts.R"))

OUT_DIR <- here("02_Proteins", "01_Differential", "c_data")
MIN_SET <- 5L

fit <- readRDS(file.path(OUT_DIR, "01_limma_DAList.rds"))
dep <- read_csv(
  file.path(OUT_DIR, "01_dep_results.csv"),
  show_col_types = FALSE
)

meta <- as.data.frame(fit$metadata)
design <- fit$design$design_matrix
contrast_matrix <- fit$design$contrast_matrix

# fry rotates the residual space and cannot take an NA, while the DEP matrix is
# the non-imputed one and is 12% missing. Set-level tests here run on the
# missForest arm, as V1's module fry and YvO's barcode panel both do. The sets
# still come from the non-imputed fit, so imputation enters the ranking only.
imputed <- readRDS(here(
  "01_Preprocess", "03_Imputation", "c_data", "DAList_imputed_missforest.rds"
))
abundance <- as.matrix(imputed$data)
stopifnot(
  identical(colnames(abundance), rownames(design)),
  identical(rownames(abundance), rownames(fit$data)),
  !anyNA(abundance)
)
correlation <- limma::duplicateCorrelation(
  abundance, design,
  block = meta$subject
)$consensus

# Mirrored pairs: each arm's signature tested against the other arm's ranking.
PAIRINGS <- tibble::tribble(
  ~set_from, ~ranked_on,
  "Training_HR", "Training_LR",
  "Training_LR", "Training_HR",
  "Acute_HR", "Acute_LR",
  "Acute_LR", "Acute_HR"
)

# Two selection criteria, as V1 reports both. Pi weights effect size and BH
# does not, so a set that holds under both is not an artefact of either.
member_rows <- function(contrast_name, criterion, direction) {
  d <- dep |> filter(.data$contrast == contrast_name)
  hit <- switch(criterion,
    pi = d$sig_pi == (if (direction == "up") 1L else -1L),
    bh = !is.na(d$adj.P.Val) & d$adj.P.Val < 0.05 &
      (if (direction == "up") d$logFC > 0 else d$logFC < 0),
    nominal = !is.na(d$P.Value) & d$P.Value < 0.05 &
      (if (direction == "up") d$logFC > 0 else d$logFC < 0)
  )
  hit[is.na(hit)] <- FALSE
  which(rownames(abundance) %in% d$uniprot_id[hit])
}

results <- purrr::pmap_dfr(PAIRINGS, function(set_from, ranked_on) {
  sets <- list()
  for (crit in c("pi", "bh", "nominal")) {
    for (dir in c("up", "down")) {
      rows <- member_rows(set_from, crit, dir)
      if (length(rows) >= MIN_SET) {
        sets[[paste(crit, dir, sep = "_")]] <- rows
      }
    }
  }
  if (!length(sets)) {
    return(tibble())
  }
  fr <- limma::fry(
    abundance,
    index = sets, design = design,
    contrast = contrast_matrix[, ranked_on],
    block = meta$subject, correlation = correlation
  )
  fr |>
    tibble::as_tibble(rownames = "set") |>
    mutate(
      set_from = set_from, ranked_on = ranked_on,
      n_proteins = lengths(sets)[.data$set], .before = 1
    )
})

summary_tbl <- results |>
  dplyr::select(
    set_from, ranked_on, set, n_proteins,
    direction = Direction, p = PValue, fdr = FDR
  ) |>
  arrange(.data$p)

write.xlsx(
  list(fry_concordance = summary_tbl),
  file.path(OUT_DIR, "03_fry_concordance.xlsx")
)

print(as.data.frame(summary_tbl), row.names = FALSE, digits = 3)
message(
  "\n", sum(summary_tbl$fdr < 0.05), " of ", nrow(summary_tbl),
  " set-ranking tests concordant at FDR < 0.05"
)

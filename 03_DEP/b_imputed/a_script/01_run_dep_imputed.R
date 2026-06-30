#!/usr/bin/env Rscript
# HRvLR COMPARISON DEP on the imputed matrices (exploratory; a_non_imputed stays the
# reported analysis). missforest = canonical (figures/WGCNA read it); imp4p and
# mscoreutils are comparison arms.
# CAVEAT: impute-before-test can inflate false positives -> logFC is checked for
# concordance against the non-imputed fit.

pacman::p_load(proteoDA, here, readr, dplyr, tibble, purrr)
source(here("03_DEP", "contrasts.R"))

CANONICAL <- "missforest"
clear_dir <- function(d) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  unlink(setdiff(list.files(d, full.names = TRUE), file.path(d, ".gitkeep")), recursive = TRUE)
}
clear_dir(here("03_DEP", "b_imputed", "c_data"))
clear_dir(here("03_DEP", "b_imputed", "b_reports"))

nd <- here("02_Normalization", "imputation", "c_data")
methods <- c(
  missforest = "DAList_imputed_missforest.rds",
  imp4p = "DAList_imputed_imp4p.rds",
  mscoreutils = "DAList_imputed_mscoreutils.rds"
)

# DEP per imputed matrix
# Same limma + duplicateCorrelation workflow as the primary arm, once per matrix.

runs <- imap(methods, function(rds, m) {
  dal <- readRDS(file.path(nd, rds))
  cat(sprintf(
    "\n[%s] imputed DEP: %d x %d%s\n", m, nrow(dal$data), ncol(dal$data),
    if (m == CANONICAL) " (canonical)" else ""
  ))
  out_dir <- here("03_DEP", "b_imputed", "c_data", m)

  dal$metadata$group <- factor(
    dal$metadata$Group_Time,
    levels = c("HR_T1", "HR_T2", "HR_T3", "LR_T1", "LR_T2", "LR_T3")
  )
  dal$metadata$subject <- dal$metadata$Subject_ID
  dal <- add_design(dal, "~ 0 + group + (1 | subject)")
  colnames(dal$design$design_matrix) <- gsub(
    "^group", "",
    colnames(dal$design$design_matrix)
  )
  dal <- add_contrasts(dal, contrasts_vector = HRVLR_CONTRASTS)
  dal <- fit_limma_model(dal)
  dal <- extract_DA_results(dal, pval_thresh = 0.10, lfc_thresh = 0, adj_method = "BH")

  ann <- as_tibble(dal$annotation) |>
    select(any_of(c("uniprot_id", "gene", "protein", "description")))
  res <- imap(dal$results, function(r, cname) {
    as_tibble(r, rownames = "uniprot_id") |>
      mutate(
        pi_score = P.Value^abs(logFC),
        sig_pi = case_when(
          pi_score < PI_THRESH & logFC > 0 ~ 1L,
          pi_score < PI_THRESH & logFC < 0 ~ -1L, TRUE ~ 0L
        ),
        contrast = cname
      ) |>
      left_join(ann, by = "uniprot_id")
  })
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write_csv(bind_rows(res), file.path(out_dir, "combined_results_pi.csv"))
  res
})

# logFC concordance vs non-imputed
# Spearman per contrast against the primary per-contrast results; high rho means
# imputation did not distort the effect estimates.

ni_dir <- here("03_DEP", "a_non_imputed", "c_data", "04_per_contrast_results")
if (dir.exists(ni_dir)) {
  ni <- map_dfr(list.files(ni_dir, pattern = "\\.csv$", full.names = TRUE), function(f) {
    read_csv(f, show_col_types = FALSE) |>
      transmute(uniprot_id,
        contrast = tools::file_path_sans_ext(basename(f)),
        logFC_ni = logFC
      )
  })
  cmp <- imap_dfr(runs, function(res, m) {
    bind_rows(res) |>
      select(uniprot_id, contrast, logFC) |>
      inner_join(ni, by = c("uniprot_id", "contrast")) |>
      group_by(contrast) |>
      summarise(
        method = m,
        rho = cor(logFC, logFC_ni, method = "spearman", use = "complete.obs"),
        n = dplyr::n(), .groups = "drop"
      )
  })
  write_csv(cmp, here("03_DEP", "b_imputed", "c_data", "logfc_concordance.csv"))
  cat("\nlogFC concordance vs non-imputed (Spearman rho):\n")
  print(cmp |> group_by(method) |> summarise(median_rho = median(rho, na.rm = TRUE)))
}

cat(sprintf(
  "\nDone: imputed DEP for %d methods -> %s\n",
  length(methods), here("03_DEP", "b_imputed", "c_data")
))

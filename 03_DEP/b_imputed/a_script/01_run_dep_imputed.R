#!/usr/bin/env Rscript
# HRvLR COMPARISON DEP on the imputed matrices (exploratory; a_non_imputed stays the
# reported analysis). missforest = canonical (figures/WGCNA read it).
#
# Read on BH, the null holds under every imputer that does not impute inside the tested factor:
# missForest (MAR), MsCoreUtils (hybrid) and Perseus (MNAR) each return zero BH<0.10 hits in all
# five HR-vs-LR and interaction contrasts, against the non-imputed arm's zero. Three contradictory
# missingness assumptions, same answer.
#
# imp4p returns 117 there, and that is circularity rather than fragility: impute.mle fits a
# separate EM inside each Group_Time cell and these contrasts test among those same cells.
# 02_imp4p_circularity.R proves it by re-imputing within random cells of equal size, which
# collapses the count back to zero. The arm is kept BECAUSE it demonstrates that.

pacman::p_load(proteoDA, here, readr, dplyr, tidyr, tibble, purrr)
source(here("03_DEP", "contrasts.R"))
source(here("04_Figures", "functions", "shared_utils.R"))

CANONICAL <- "missforest"
clear_dir(here("03_DEP", "b_imputed", "c_data"))
clear_dir(here("03_DEP", "b_imputed", "b_reports"))

nd <- here("02_Normalization", "imputation", "c_data")
methods <- c(
  missforest = "DAList_imputed_missforest.rds",
  imp4p = "DAList_imputed_imp4p.rds",
  mscoreutils = "DAList_imputed_mscoreutils.rds",
  perseus = "DAList_imputed_perseus.rds"
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
  dal <- add_contrasts(dal, contrasts_vector = HRVLR_CONTRASTS)
  dal <- fit_limma_model(dal)
  dal <- extract_DA_results(dal, pval_thresh = 0.10, lfc_thresh = 0, adj_method = "BH")

  ann <- as_tibble(dal$annotation) |>
    select(any_of(c("uniprot_id", "gene", "protein", "description")))
  res <- imap(dal$results, function(r, cname) {
    as_tibble(r, rownames = "uniprot_id") |>
      add_pi_score() |>
      mutate(contrast = cname) |>
      left_join(ann, by = "uniprot_id")
  })
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write_csv(bind_rows(res), file.path(out_dir, "combined_results_pi.csv"))
  res
})

# Sensitivity table: BH beside Pi, per contrast per arm.
# BH is the metric the headline rests on. Pi is reported alongside it but must never carry a
# cross-arm robustness claim: Pi = p^|log2FC| exponentiates the fold change, which is precisely
# what imputation distorts. Perseus is the proof — the MOST Pi-hits of any arm and the SAME BH
# count as the non-imputed fit.

NULL_CONTRASTS <- c(
  "Baseline_HRvLR", "Trained_HRvLR", "Acute_HRvLR",
  "Training_Interaction", "Acute_Interaction"
)

ni_dir <- here("03_DEP", "a_non_imputed", "c_data", "04_per_contrast_results")
ni <- map_dfr(list.files(ni_dir, pattern = "\\.csv$", full.names = TRUE), function(f) {
  read_csv(f, show_col_types = FALSE) |>
    add_pi_score() |>
    transmute(uniprot_id,
      contrast = tools::file_path_sans_ext(basename(f)),
      logFC_ni = logFC, adj_ni = adj.P.Val, pi_ni = pi_score
    )
})

arm_counts <- function(d, m) {
  d |>
    group_by(contrast) |>
    summarise(
      method = m,
      n_BH_10 = sum(adj.P.Val < 0.10, na.rm = TRUE),
      n_BH_05 = sum(adj.P.Val < 0.05, na.rm = TRUE),
      n_pi = sum(pi_score < PI_THRESH, na.rm = TRUE),
      .groups = "drop"
    )
}

sens <- bind_rows(
  arm_counts(
    ni |> transmute(contrast, adj.P.Val = adj_ni, pi_score = pi_ni),
    "non_imputed"
  ),
  imap_dfr(runs, \(res, m) arm_counts(bind_rows(res), m))
) |>
  mutate(is_null_contrast = contrast %in% NULL_CONTRASTS) |>
  arrange(method, contrast)

write_csv(sens, here("03_DEP", "b_imputed", "c_data", "sensitivity_bh_vs_pi.csv"))

cat("\nSensitivity — BH beside Pi, summed over the five HR-vs-LR / interaction contrasts:\n")
sens |>
  filter(is_null_contrast) |>
  group_by(method) |>
  summarise(BH_10 = sum(n_BH_10), BH_05 = sum(n_BH_05), Pi = sum(n_pi), .groups = "drop") |>
  arrange(BH_10) |>
  as.data.frame() |>
  print()

# logFC concordance vs non-imputed, STRATIFIED BY MISSINGNESS.
# The unstratified rho is worthless here: it is dominated by the ~944 proteins with no missing
# values, where every imputer is a no-op by construction. It ranked imp4p 2nd best (0.928) while
# imp4p was destroying the null, and Perseus worst (0.491) while Perseus preserved it exactly.
# Imputation can only move a protein that had a gap, so the rho must be read within gap strata.

n_na <- rowSums(is.na(readRDS(here("02_Normalization", "c_data", "DAList_normalized.rds"))$data))
miss_strata <- tibble(
  uniprot_id = names(n_na),
  n_missing = as.integer(n_na),
  stratum = cut(n_na,
    breaks = c(-1, 0, 9, Inf),
    labels = c("0 (imputation is a no-op)", "1-9", ">=10")
  )
)

cmp <- imap_dfr(runs, function(res, m) {
  bind_rows(res) |>
    select(uniprot_id, contrast, logFC) |>
    inner_join(ni, by = c("uniprot_id", "contrast")) |>
    inner_join(miss_strata, by = "uniprot_id") |>
    group_by(contrast, stratum) |>
    summarise(
      method = m,
      rho = cor(logFC, logFC_ni, method = "spearman", use = "complete.obs"),
      n = dplyr::n(), .groups = "drop"
    )
})
write_csv(cmp, here("03_DEP", "b_imputed", "c_data", "logfc_concordance.csv"))

cat("\nlogFC concordance vs non-imputed (median Spearman rho by missingness stratum):\n")
cmp |>
  group_by(method, stratum) |>
  summarise(median_rho = round(median(rho, na.rm = TRUE), 3), .groups = "drop") |>
  pivot_wider(names_from = stratum, values_from = median_rho) |>
  as.data.frame() |>
  print()

cat(sprintf(
  "\nDone: imputed DEP for %d methods -> %s\n",
  length(methods), here("03_DEP", "b_imputed", "c_data")
))

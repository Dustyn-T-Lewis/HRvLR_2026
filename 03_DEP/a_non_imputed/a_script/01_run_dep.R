# HRvLR Differential Expression — limma pipeline (proteoDA)
# Design: 2x3 factorial (Responder x Timepoint), repeated measures on subject
# T1 = baseline, T2 = 72hr post-training, T3 = 1hr acute post-bout
# Input:  cycloess-normalized, non-imputed (limma handles NAs per-protein)
#
# Contrasts: nine, defined in 03_DEP/contrasts.R (HRVLR_CONTRASTS)
#
# References:
#   Ritchie et al. 2015, Nucleic Acids Res 43(7):e47 — limma
#   Smyth 2005, Stat Appl Genet Mol Biol 3(1):3 — duplicateCorrelation
#   Karpievitch et al. 2012, BMC Bioinform 13(S16):S5 — non-imputed limma
#   Xiao et al. 2014, Bioinformatics 30(6):801-807 — Pi-score
#     Pi = p^|logFC|; threshold Pi < 0.05 <-> original pi > 1.3

pacman::p_load(dplyr, tidyr, tibble, readr, purrr, proteoDA, here)
source(here("03_DEP", "contrasts.R"))

cfg <- list(
  norm_csv = here("02_Normalization", "c_data", "normalized.csv"),
  norm_rds = here("02_Normalization", "c_data", "DAList_normalized.rds"),
  data_dir = here("03_DEP", "a_non_imputed", "c_data"),
  report_dir = here("03_DEP", "a_non_imputed", "b_reports"),
  proteoDA_dir = here("03_DEP", "a_non_imputed", "b_reports", "01_proteoDA"),
  pval_thresh = 0.10,
  lfc_thresh = 0,
  adj_method = "BH",
  pi_thresh = PI_THRESH
)

dir.create(cfg$data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$proteoDA_dir, recursive = TRUE, showWarnings = FALSE)

required_meta_cols <- c("Col_ID", "Subject_ID", "Group", "Timepoint", "Group_Time")

# load data and build metadata

df <- read_csv(cfg$norm_csv, show_col_types = FALSE)

ann_cols <- c("uniprot_id", "protein", "gene", "description")
ann <- df[, ann_cols]
samp_names <- setdiff(names(df), ann_cols)
mat <- as.matrix(df[, samp_names])
rownames(mat) <- ann$uniprot_id

cat(sprintf(
  "Loaded: %d proteins x %d samples | missing: %d (%.1f%%)\n",
  nrow(mat), ncol(mat), sum(is.na(mat)),
  100 * sum(is.na(mat)) / length(mat)
))

# Canonical metadata from normalisation DAList (not regex-derived)
dal_norm <- readRDS(cfg$norm_rds)
dal_meta <- as.data.frame(dal_norm$metadata)
missing_meta_cols <- setdiff(required_meta_cols, names(dal_meta))
if (length(missing_meta_cols)) {
  stop(sprintf(
    "Missing required metadata columns in normalized DAList: %s",
    paste(missing_meta_cols, collapse = ", ")
  ))
}

meta <- tibble(
  sample_id  = dal_meta$Col_ID,
  responder  = dal_meta$Group,
  time       = dal_meta$Timepoint,
  group      = dal_meta$Group_Time,
  subject    = dal_meta$Subject_ID
)
meta$responder <- factor(meta$responder, levels = c("HR", "LR"))
meta$time <- factor(meta$time, levels = c("T1", "T2", "T3"))
meta$group <- factor(meta$group,
  levels = c(
    "HR_T1", "HR_T2", "HR_T3",
    "LR_T1", "LR_T2", "LR_T3"
  )
)

if (any(is.na(meta$subject)) || any(meta$subject == "")) {
  stop("Subject_ID must be present for duplicateCorrelation blocking.")
}

print(table(meta$responder, meta$time))
stopifnot(setequal(colnames(mat), meta$sample_id))

# build the DAList

meta_df <- as.data.frame(meta)
rownames(meta_df) <- meta$sample_id

dal <- DAList(
  data       = mat,
  annotation = as.data.frame(ann),
  metadata   = meta_df,
  tags       = list(norm_method = "cycloess", normalized = TRUE)
)

# design

dal <- add_design(dal, "~ 0 + group + (1 | subject)")

# contrasts

dal <- add_contrasts(dal, contrasts_vector = HRVLR_CONTRASTS)

# fit and extract

dal <- fit_limma_model(dal)

within_cor <- dal$eBayes_fit$correlation %||%
  dal$tags$duplicate_correlation %||% NA_real_
if (!is.na(within_cor)) cat(sprintf("Within-subject correlation: %.3f\n", within_cor))

# FDR 0.10: appropriate for exploratory proteomics with small n
# (Benjamini & Hochberg 1995 used 0.10 in examples; Choi et al. 2008,
# BMC Bioinform 9:43, recommend FDR 0-15% for label-free proteomics).
# Pi-score (Pi < 0.05) provides a secondary effect-size-weighted filter.
dal <- extract_DA_results(dal,
  pval_thresh = cfg$pval_thresh,
  lfc_thresh  = cfg$lfc_thresh,
  adj_method  = cfg$adj_method
)

# add the pi-score to the results
# Pi = p^|logFC| (Xiao et al. 2014); lower = more significant
contrast_names <- names(dal$results)

for (cname in contrast_names) {
  res <- dal$results[[cname]]
  res$pi_score <- res$P.Value^abs(res$logFC)
  res$sig_pi <- case_when(
    res$pi_score < cfg$pi_thresh & res$logFC > 0 ~ 1L,
    res$pi_score < cfg$pi_thresh & res$logFC < 0 ~ -1L,
    TRUE ~ 0L
  )
  dal$results[[cname]] <- res
}

# proteoDA Excel formatting expects gene_symbol column
dal$annotation$gene_symbol <- dal$annotation$gene

# save the fitted DAList

saveRDS(dal, file.path(cfg$data_dir, "01_limma_DAList.rds"))

# proteoDA interactive report (non-essential; warn on failure)

tryCatch(
  write_limma_plots(dal,
    grouping_column = "group",
    output_dir      = cfg$proteoDA_dir,
    table_columns   = c("uniprot_id", "gene", "protein"),
    title_column    = "gene",
    overwrite       = TRUE
  ),
  error = function(e) warning("write_limma_plots failed: ", conditionMessage(e), call. = FALSE)
)

# write tables
# Per-contrast CSVs, combined results (wide), and formatted Excel workbook
# with conditional formatting & UniProt hyperlinks — all handled by proteoDA.
# Pi-score columns pass through because they were injected into dal$results above.

write_limma_tables(dal,
  output_dir        = cfg$data_dir,
  contrasts_subdir  = "04_per_contrast_results",
  summary_csv       = "02_DA_summary_base.csv",
  combined_file_csv = "03_combined_results.csv",
  spreadsheet_xlsx  = "05_results.xlsx",
  annot_cols        = c("uniprot_id", "gene", "protein", "description"),
  overwrite         = TRUE
)

# proteoDA summary only has sig.PVal/sig.FDR; replace with Pi-enriched version
file.remove(file.path(cfg$data_dir, "02_DA_summary_base.csv"))

# DA summary
# proteoDA's summarize_contrast_DA lacks Pi-score columns, so we build our own.

da_count_row <- function(res, cname, type) {
  if (type == "nonsig") {
    n_pval <- sum(res$P.Value >= cfg$pval_thresh, na.rm = TRUE)
    n_fdr <- sum(res$adj.P.Val >= cfg$pval_thresh, na.rm = TRUE)
    n_pi <- sum(res$sig_pi == 0, na.rm = TRUE)
    n_05 <- sum(res$adj.P.Val >= 0.05, na.rm = TRUE)
    n_10 <- sum(res$adj.P.Val >= 0.10, na.rm = TRUE)
  } else {
    dir <- if (type == "up") res$logFC > 0 else res$logFC < 0
    pi_target <- if (type == "up") 1L else -1L
    n_pval <- sum(res$P.Value < cfg$pval_thresh & dir, na.rm = TRUE)
    n_fdr <- sum(res$adj.P.Val < cfg$pval_thresh & dir, na.rm = TRUE)
    n_pi <- sum(res$sig_pi == pi_target, na.rm = TRUE)
    n_05 <- sum(res$adj.P.Val < 0.05 & dir, na.rm = TRUE)
    n_10 <- sum(res$adj.P.Val < 0.10 & dir, na.rm = TRUE)
  }
  tibble(
    contrast = cname, type = type,
    sig.PVal = n_pval, sig.FDR = n_fdr,
    pval_thresh = cfg$pval_thresh, lfc_thresh = cfg$lfc_thresh,
    p_adj_method = cfg$adj_method,
    sig.Pi = n_pi, sig.FDR.05 = n_05, sig.FDR.10 = n_10
  )
}

da_summary <- map_dfr(contrast_names, function(cname) {
  res <- dal$results[[cname]]
  map_dfr(c("up", "down", "nonsig"), function(type) da_count_row(res, cname, type))
})

write_csv(da_summary, file.path(cfg$data_dir, "02_DA_summary.csv"))

# summary

print(dal$design$contrast_matrix)
print(da_summary)

cat(sprintf(
  "Done: 01_run_dep.R — %d contrasts -> %s/\n",
  length(contrast_names), cfg$data_dir
))

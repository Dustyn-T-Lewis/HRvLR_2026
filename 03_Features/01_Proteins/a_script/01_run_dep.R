# HRvLR Differential Expression — limma pipeline (proteoDA)
# Input: cycloess-normalized, non-imputed (limma handles NAs per-protein)
# T1 = baseline, T2 = 72hr post-training, T3 = 1hr acute post-bout
#
# Two contrast families, one script, one matrix, one estimator:
#
#   categorical  ~ 0 + group + (1 | subject)      the nine HRVLR_CONTRASTS
#                2x3 factorial, Responder x Timepoint
#   pooled       ~ 0 + time + (1 | subject)       the two POOLED_CONTRASTS
#                no arm term anywhere, averaged across all 16 subjects
#
# Both are defined in 03_Features/contrasts.R and fitted here through proteoDA
# so the committed protein-level numbers and a feature_design() refit can be
# checked against each other (verify_protein_equivalence()). They ran as two
# directories until 2026-08-31; the designs differ, so proteoDA necessarily
# builds a DAList apiece, but the combined results table they feed is one file
# carrying every contrast, which is what stops the two drifting apart.
#
# References:
#   Ritchie et al. 2015, Nucleic Acids Res 43(7):e47 — limma
#   Smyth, Michaud & Scott 2005, Bioinformatics 21(9):2067 —
#   duplicateCorrelation
#   Smyth 2004, Stat Appl Genet Mol Biol 3:1 — empirical Bayes moderation (the
#   eBayes engine)
#   Phipson et al. 2016, Ann Appl Stat 10(2):946 — robust empirical Bayes
#   Xiao et al. 2014, Bioinformatics 30(6):801-807 — Pi-score
#     Pi = p^|logFC|; threshold Pi < 0.05 <-> original pi > 1.3. A transformed
#     raw p, bounded in [0,1], controlling no error rate. See
#     04_pi_permutation.R.
#
# On running limma without imputing: Karpievitch et al. 2012, BMC Bioinform
# 13(S16):S5 is the source of the known COST of this choice, not a licence for
# it — it warns that complete-case analysis yields downward-biased standard
# errors, i.e. it is anti-conservative. We accept that and report it, because a
# method biased toward false positives returning zero BH hits in every HR-vs-LR
# contrast makes the null stronger, not weaker. The imputed arms in
# 03_Features/01_Proteins/imputed are the robustness check.

pacman::p_load(dplyr, tibble, readr, purrr, proteoDA, here)
source(here("03_Features", "contrasts.R"))

cfg <- list(
  norm_csv = here("02_Normalization", "c_data", "normalized.csv"),
  norm_rds = here("02_Normalization", "c_data", "DAList_normalized.rds"),
  data_dir = here("03_Features", "01_Proteins", "c_data"),
  report_dir = here("03_Features", "01_Proteins", "b_reports"),
  proteoDA_dir = here(
    "03_Features", "01_Proteins", "b_reports", "01_proteoDA"
  ),
  pval_thresh = 0.10,
  lfc_thresh = 0,
  adj_method = "BH",
  pi_thresh = PI_THRESH
)

# Each family names the design column it groups on, the formula, and its
# contrast vector. Everything downstream of the DAList is identical, so the
# fit runs once per entry and nothing else forks.
FAMILIES <- list(
  categorical = list(
    grouping = "group",
    formula = "~ 0 + group + (1 | subject)",
    contrasts = HRVLR_CONTRASTS
  ),
  pooled = list(
    grouping = "time",
    formula = "~ 0 + time + (1 | subject)",
    contrasts = POOLED_CONTRASTS
  )
)

dir.create(cfg$data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$proteoDA_dir, recursive = TRUE, showWarnings = FALSE)

required_meta_cols <- c(
  "Col_ID", "Subject_ID", "Group", "Timepoint", "Group_Time"
)

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
meta$group <- factor(meta$group, levels = GROUP_LEVELS)

if (any(is.na(meta$subject)) || any(meta$subject == "")) {
  stop("Subject_ID must be present for duplicateCorrelation blocking.")
}

print(table(meta$responder, meta$time))
stopifnot(setequal(colnames(mat), meta$sample_id))

meta_df <- as.data.frame(meta)
rownames(meta_df) <- meta$sample_id

# Selection is by BH. Pi and raw p are reported beside it; neither selects.
# BH at 0.10 is our threshold for exploratory n=16 proteomics, not a
# literature-mandated one. BH is applied WITHIN each contrast (topTable is
# called per coef, and decideTests defaults to method = "separate"), never
# across the nine, and never across the two families.
#
# Pi used to be the selection criterion. 04_pi_permutation.R retired it:
# shuffling the arm label across subjects produces MORE Pi hits than the real
# labels do (235 observed against a permuted median of 274 across the five
# HR-vs-LR and interaction contrasts, no contrast below emp p = 0.22). Never
# quote a sig.Pi count without that null beside it. That permutation was run
# on the arm label; whether it holds for a shuffled timepoint label in the
# pooled family is untested, so read pooled Pi counts with that gap stated.
da_count_row <- function(res, cname, type) {
  if (type == "nonsig") {
    n_p10 <- sum(res$P.Value >= 0.10, na.rm = TRUE)
    n_pi <- sum(res$sig_pi == 0, na.rm = TRUE)
    n_q05 <- sum(res$adj.P.Val >= 0.05, na.rm = TRUE)
    n_q10 <- sum(res$adj.P.Val >= 0.10, na.rm = TRUE)
  } else {
    dir <- if (type == "up") res$logFC > 0 else res$logFC < 0
    pi_target <- if (type == "up") 1L else -1L
    n_p10 <- sum(res$P.Value < 0.10 & dir, na.rm = TRUE)
    n_pi <- sum(res$sig_pi == pi_target, na.rm = TRUE)
    n_q05 <- sum(res$adj.P.Val < 0.05 & dir, na.rm = TRUE)
    n_q10 <- sum(res$adj.P.Val < 0.10 & dir, na.rm = TRUE)
  }
  tibble(
    contrast = cname, type = type,
    sig.P.10 = n_p10, sig.FDR.05 = n_q05, sig.FDR.10 = n_q10, sig.Pi = n_pi,
    lfc_thresh = cfg$lfc_thresh, p_adj_method = cfg$adj_method
  )
}

run_family <- function(family) {
  spec <- FAMILIES[[family]]
  cat(sprintf(
    "\n=== %s family: %d contrasts ===\n",
    family, length(spec$contrasts)
  ))

  dal <- DAList(
    data       = mat,
    annotation = as.data.frame(ann),
    metadata   = meta_df,
    tags       = list(norm_method = "cycloess", normalized = TRUE)
  )
  dal <- add_design(dal, spec$formula)
  dal <- add_contrasts(dal, contrasts_vector = spec$contrasts)
  dal <- fit_limma_model(dal)

  within_cor <- limma::duplicateCorrelation(
    dal$data, dal$design$design_matrix,
    block = meta$subject
  )$consensus
  cat(sprintf("Within-subject correlation: %.3f\n", within_cor))

  dal <- extract_DA_results(dal,
    pval_thresh = cfg$pval_thresh,
    lfc_thresh  = cfg$lfc_thresh,
    adj_method  = cfg$adj_method
  )

  # add_pi_score is defined beside the threshold in 03_Features/contrasts.R so
  # both DEP arms gate identically.
  contrast_names <- names(dal$results)
  for (cname in contrast_names) {
    dal$results[[cname]] <- add_pi_score(dal$results[[cname]], cfg$pi_thresh)
  }

  # proteoDA Excel formatting expects gene_symbol column
  dal$annotation$gene_symbol <- dal$annotation$gene

  saveRDS(dal, file.path(
    cfg$data_dir, sprintf("01_limma_DAList_%s.rds", family)
  ))

  # proteoDA interactive report (non-essential; warn on failure)
  tryCatch(
    write_limma_plots(dal,
      grouping_column = spec$grouping,
      output_dir      = file.path(cfg$proteoDA_dir, family),
      table_columns   = c("uniprot_id", "gene", "protein"),
      title_column    = "gene",
      overwrite       = TRUE
    ),
    error = function(e) {
      warning("write_limma_plots failed: ", conditionMessage(e), call. = FALSE)
    }
  )

  # Per-contrast CSVs, combined results (wide) and the formatted workbook, all
  # handled by proteoDA. Only the wide combined CSV carries the pi columns
  # injected into dal$results above; the per-contrast CSVs and the workbook
  # sheets are built from a fixed column list and stop at sig.FDR. Contrast
  # names never collide across families, so the per-contrast CSVs share one
  # directory.
  write_limma_tables(dal,
    output_dir        = cfg$data_dir,
    contrasts_subdir  = "04_per_contrast_results",
    summary_csv       = sprintf("base_summary_%s.csv", family),
    combined_file_csv = sprintf("combined_%s.csv", family),
    spreadsheet_xlsx  = sprintf("05_results_%s.xlsx", family),
    annot_cols        = c("uniprot_id", "gene", "protein", "description"),
    overwrite         = TRUE
  )
  # proteoDA's summary has only sig.PVal/sig.FDR and drives both off the single
  # pval_thresh argument, so its `sig.PVal` means p < 0.10 and its `sig.FDR`
  # duplicates the q < 0.10 column. Ours below names every threshold it applies.
  file.remove(file.path(cfg$data_dir, sprintf("base_summary_%s.csv", family)))

  summary_tbl <- map_dfr(contrast_names, function(cname) {
    res <- dal$results[[cname]]
    map_dfr(
      c("up", "down", "nonsig"),
      function(type) da_count_row(res, cname, type)
    )
  }) |>
    mutate(family = family, .before = "contrast")

  print(dal$design$contrast_matrix)
  list(summary = summary_tbl, n_contrasts = length(contrast_names))
}

fits <- lapply(setNames(names(FAMILIES), names(FAMILIES)), run_family)

# One combined table for every downstream reader. proteoDA writes one wide CSV
# per family (rows = proteins, columns = statistic x contrast). Annotation and
# per-sample intensity columns are identical across the two -- same matrix, same
# row order -- so the merge asserts that and then appends only the columns each
# later family adds. A join on annotation alone would duplicate all 45 intensity
# columns as .x/.y.
tables <- lapply(names(FAMILIES), function(family) {
  read_csv(
    file.path(cfg$data_dir, sprintf("combined_%s.csv", family)),
    show_col_types = FALSE
  )
})
stopifnot(all(vapply(
  tables[-1], function(t) identical(t$uniprot_id, tables[[1]]$uniprot_id),
  logical(1)
)))
combined <- Reduce(
  function(a, b) bind_cols(a, b[setdiff(names(b), names(a))]),
  tables
)
write_csv(combined, file.path(cfg$data_dir, "03_combined_results.csv"))
file.remove(file.path(
  cfg$data_dir, sprintf("combined_%s.csv", names(FAMILIES))
))

da_summary <- bind_rows(lapply(fits, `[[`, "summary"))
write_csv(da_summary, file.path(cfg$data_dir, "02_DA_summary.csv"))
print(as.data.frame(da_summary))

cat(sprintf(
  "Done: 01_run_dep.R — %d contrasts across %d families -> %s/\n",
  sum(vapply(fits, `[[`, integer(1), "n_contrasts")), length(FAMILIES),
  cfg$data_dir
))

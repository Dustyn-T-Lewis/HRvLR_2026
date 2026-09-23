# The nine HR-versus-LR contrasts, fitted through proteoDA.
#
# Four within-arm changes, three between-arm differences at each timepoint, and
# two interactions. The design is a means model over the six Group_Time cells
# with subject as the blocking factor, so the repeated measures are carried by
# duplicateCorrelation rather than ignored.
#
# This reproduces V1's committed numbers exactly; 02_verify_v1.R is the proof.
# The point of refitting rather than copying is that everything downstream in
# this project reads one DAList, and a stage that cannot rebuild its own input
# is a stage nobody can check.

pacman::p_load(
  here, dplyr, tidyr, tibble, readr, purrr, limma, proteoDA, openxlsx
)

source(here("functions", "contrasts.R"))
source(here("functions", "association.R"))

OUT_DIR <- here("02_Proteins", "01_Differential", "c_data")
RPT_DIR <- here("02_Proteins", "01_Differential", "b_reports")

dal <- readRDS(here(
  "01_Preprocess", "02_Normalization", "c_data", "DAList_normalized.rds"
))
src <- as.data.frame(dal$metadata)
stopifnot(identical(src$Col_ID, colnames(dal$data)))

# proteoDA reads the design formula's terms from metadata column names, and
# carries the random effect as (1 | subject) inside that formula rather than as
# an argument. The renaming below is what makes the formula readable.
meta <- data.frame(
  sample_id = src$Col_ID,
  responder = factor(src$Group, levels = c("HR", "LR")),
  time = factor(src$Timepoint, levels = c("T1", "T2", "T3")),
  group = factor(src$Group_Time, levels = GROUP_LEVELS),
  subject = src$Subject_ID
)
rownames(meta) <- meta$sample_id
stopifnot(!any(is.na(meta$subject) | meta$subject == ""))
dal$metadata <- meta

dal <- dal |>
  proteoDA::add_design("~ 0 + group + (1 | subject)") |>
  proteoDA::add_contrasts(contrasts_vector = HRVLR_CONTRASTS)

fit <- proteoDA::fit_limma_model(dal)

# proteoDA computes the consensus correlation internally and discards it, so it
# is recomputed here for the record. A value near zero would mean the blocking
# is doing nothing; here it is 0.189.
within_cor <- limma::duplicateCorrelation(
  dal$data, fit$design$design_matrix,
  block = meta$subject
)$consensus

res <- proteoDA::extract_DA_results(
  fit,
  pval_thresh = 0.05, lfc_thresh = 0, adj_method = "BH"
)

# extract_DA_results() returns statistics keyed only by rownames. Those
# rownames are the annotation's rownames, so the two align by position; the
# stopifnot is what makes relying on that safe rather than hopeful. Carrying
# the ids as columns lets every downstream stage join instead of trusting order.
annotation <- fit$annotation[
  , c("uniprot_id", "gene", "protein", "description")
]
stopifnot(identical(rownames(res$results[[1]]), rownames(annotation)))

combined <- res$results |>
  lapply(function(d) dplyr::bind_cols(annotation, d)) |>
  dplyr::bind_rows(.id = "contrast") |>
  tibble::as_tibble() |>
  add_pi_score()

# Three significance layers per contrast, reported side by side as V1 did.
# The pi-score is the primary criterion (Xiao Eq. 2, pi = p^|log2FC|), BH is
# the conservative one, and the nominal count is the exploratory one. Proteins
# with an inestimable cell mean come back NA and are never tested, so the BH
# denominator is below 1900 and differs by contrast; recording it stops anyone
# reading a q as if it spanned everything.
contrast_summary <- combined |>
  summarise(
    n_tested = sum(!is.na(.data$P.Value)),
    n_untested = sum(is.na(.data$P.Value)),
    n_nominal = sum(.data$P.Value < 0.05, na.rm = TRUE),
    n_pi = sum(.data$sig_pi != 0L, na.rm = TRUE),
    n_pi_up = sum(.data$sig_pi == 1L, na.rm = TRUE),
    n_pi_down = sum(.data$sig_pi == -1L, na.rm = TRUE),
    n_bh05 = sum(.data$adj.P.Val < 0.05, na.rm = TRUE),
    min_bh = min(.data$adj.P.Val, na.rm = TRUE),
    .by = "contrast"
  ) |>
  arrange(match(.data$contrast, CONTRAST_NAMES))

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_DIR, recursive = TRUE, showWarnings = FALSE)

saveRDS(fit, file.path(OUT_DIR, "01_limma_DAList.rds"))
write_csv(combined, file.path(OUT_DIR, "01_dep_results.csv"))
write_csv(contrast_summary, file.path(RPT_DIR, "contrast_summary.csv"))
write.xlsx(
  c(
    list(contrast_summary = contrast_summary),
    split(combined, factor(combined$contrast, levels = CONTRAST_NAMES))
  ),
  file.path(OUT_DIR, "01_dep_results.xlsx")
)

proteoDA::write_limma_plots(
  res,
  grouping_column = "group", output_dir = RPT_DIR,
  table_columns = c("uniprot_id", "gene", "protein"), title_column = "gene",
  overwrite = TRUE
)

# The handoff every later protein-level step reads. Statistics come as protein
# by contrast matrices so a pathway or network stage can index them without
# re-deriving the fit; the three subject windows are the inputs the
# classification and association screens share across feature levels.
stat_matrix <- function(col) {
  combined |>
    select("uniprot_id", "contrast", value = all_of(col)) |>
    pivot_wider(names_from = "contrast", values_from = "value") |>
    tibble::column_to_rownames("uniprot_id") |>
    as.matrix()
}
abund <- as.matrix(fit$data)
sample_meta <- sample_metadata()
stopifnot(identical(sample_meta$sample_id, colnames(abund)))
proteins <- list(
  abund = abund,
  meta = sample_meta,
  annotation = annotation,
  stats = set_names(c("logFC", "t", "P.Value", "adj.P.Val")) |>
    map(stat_matrix),
  windows = set_names(c("T1", "training", "acute")) |>
    map(\(w) subject_window(abund, w, sample_meta))
)
saveRDS(proteins, file.path(OUT_DIR, "proteins.rds"))

print(as.data.frame(contrast_summary), row.names = FALSE, digits = 3)
message(
  "\nwithin-subject correlation: ", round(within_cor, 4),
  " | pi-score hits: ", sum(contrast_summary$n_pi),
  " | nominal: ", sum(contrast_summary$n_nominal),
  " | BH: ", sum(contrast_summary$n_bh05)
)

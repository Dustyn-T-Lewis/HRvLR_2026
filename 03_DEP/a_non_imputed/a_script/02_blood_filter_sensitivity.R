#!/usr/bin/env Rscript
# Does the null depend on having deleted the blood-tagged proteins?
#
# The primary arm removes 120 proteins on HPA's "Secreted to blood" tag unless they clear a
# muscle rescue (myonuclei >= 20 AND blood_conc < 1e9). That rule deletes 5% of the proteome,
# and the deleted set is enriched BY CONSTRUCTION for the secreted / ECM / myokine compartment.
# So "no HR-vs-LR difference in the muscle secretome" is currently an untested claim rather than
# a null result. This tests it.
#
# The rescue is also known to be imperfect. HPA's tag is a sequence-based prediction (signal
# peptide), not a measurement of carryover, and the transcript floor overrides the measured blood
# concentration: C1QBP (mitochondrial, blood 1.5e6 pg/L), HMGB2 (nuclear, 8.3e5), PPIB, LGALS1
# and CTSB (an exercise myokine) are all deleted despite plasma levels 4-5 orders of magnitude
# below the genuine plasma proteins.
#
# Sensitivity arm: add the 120 back, change nothing else, refit the same nine contrasts.
# The curated list stays removed - keratins, globins and the 95 blood proteins matched by
# accession in 00_input/blood_contaminants.csv, none of which is muscle under any reading.
#
# NOTE: the protein count drives the normalization. limma's adaptive.span sets the loess span
# from nrow, so a larger matrix normalizes on a different span (0.5075 -> reported below). That
# is a genuine consequence of the filter, not an artifact of this script.

pacman::p_load(proteoDA, here, readxl, readr, dplyr, stringr, tibble, purrr)
source(here("03_DEP", "contrasts.R"))
source(here("01_Filtering", "a_script", "filter_config.R"))

NULL_CONTRASTS <- c(
  "Baseline_HRvLR", "Trained_HRvLR", "Acute_HRvLR",
  "Training_Interaction", "Acute_Interaction"
)

calls <- read_csv(here("01_Filtering", "c_data", "protein_calls.csv"), show_col_types = FALSE)
primary <- readRDS(here("01_Filtering", "c_data", "DAList_filtered.rds"))
kept_samples <- colnames(primary$data)

raw <- read_excel(here("00_input", "HRvLR_raw.xlsx"))
annot_cols <- c("uniprot_id", "protein", "gene", "description", "n_seq")
metadata <- as.data.frame(read_csv(here("00_input", "HRvLR_meta.csv"), show_col_types = FALSE))
rownames(metadata) <- metadata$Col_ID

# Everything the primary arm kept, plus the proteins it removed on the blood tag alone.
readmit <- calls$uniprot_id[calls$verdict == "remove: plasma"]
keep_ids <- c(calls$uniprot_id[!str_starts(calls$verdict, "remove")], readmit)

annotation <- as.data.frame(raw[raw$uniprot_id %in% keep_ids, annot_cols])
intensity <- as.data.frame(data.matrix(raw[raw$uniprot_id %in% keep_ids, metadata$Col_ID]))
rownames(intensity) <- rownames(annotation) <- annotation$uniprot_id

dal <- zero_to_missing(DAList(
  data = intensity, annotation = annotation, metadata = metadata
))
dal <- filter_samples(dal, Col_ID %in% kept_samples)
dal <- filter_proteins_by_group(dal,
  min_reps = filter_cfg$miss_min_reps,
  min_groups = filter_cfg$miss_min_groups,
  grouping_column = "Group_Time"
)
dal <- normalize_data(dal, norm_method = "cycloess")

cat(sprintf(
  "readmitted %d blood-tagged proteins | %d x %d (primary: %d x %d) | loess span %.4f\n",
  length(readmit), nrow(dal$data), ncol(dal$data),
  nrow(primary$data), length(kept_samples),
  limma::chooseLowessSpan(nrow(dal$data))
))

dal$metadata$group <- factor(
  dal$metadata$Group_Time,
  levels = c("HR_T1", "HR_T2", "HR_T3", "LR_T1", "LR_T2", "LR_T3")
)
dal$metadata$subject <- dal$metadata$Subject_ID
dal <- add_design(dal, "~ 0 + group + (1 | subject)")
dal <- add_contrasts(dal, contrasts_vector = HRVLR_CONTRASTS)
dal <- fit_limma_model(dal)
dal <- extract_DA_results(dal, pval_thresh = 0.10, lfc_thresh = 0, adj_method = "BH")

res <- imap_dfr(dal$results, function(r, cname) {
  as_tibble(r, rownames = "uniprot_id") |>
    mutate(contrast = cname, readmitted = uniprot_id %in% readmit)
}) |>
  left_join(select(calls, uniprot_id, gene), by = "uniprot_id")

summary_tbl <- res |>
  group_by(contrast) |>
  summarise(
    n_tested = sum(!is.na(adj.P.Val)),
    BH_10 = sum(adj.P.Val < 0.10, na.rm = TRUE),
    BH_10_readmitted = sum(adj.P.Val < 0.10 & readmitted, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(is_null_contrast = contrast %in% NULL_CONTRASTS)

report_dir <- here("03_DEP", "a_non_imputed", "b_reports")
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
write_csv(summary_tbl, file.path(report_dir, "blood_filter_sensitivity.csv"))

hits <- res |>
  filter(adj.P.Val < 0.10, contrast %in% NULL_CONTRASTS) |>
  arrange(adj.P.Val) |>
  select(contrast, gene, uniprot_id, logFC, adj.P.Val, readmitted)
write_csv(hits, file.path(report_dir, "blood_filter_sensitivity_hits.csv"))

print(as.data.frame(summary_tbl), row.names = FALSE)
cat(sprintf(
  "\nBH<0.10 across the five HR-vs-LR / interaction contrasts: %d (primary arm: 0)\n",
  sum(summary_tbl$BH_10[summary_tbl$is_null_contrast])
))
cat("If this is 0, the null does not depend on the blood filter and the secretome claim is earned.\n")

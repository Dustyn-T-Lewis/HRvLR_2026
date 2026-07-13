#!/usr/bin/env Rscript
# HRvLR Stage 01: dedup -> contaminants + HPA blood removal -> outlier consensus ->
# group-wise missingness. Writes the filtered, un-normalized DAList for Stage 02, a verdict
# for every protein, and per-sample contamination indices.
#
# Removal is decided once, in protein_calls$verdict. Everything else is a view of that table.
# Blood logic follows CvH: remove iff blood-derived (secreted-to-blood | immunoglobulin |
# erythrocyte-high) AND NOT rescued as muscle. Absence from HPA is UNKNOWN, never a reason to
# remove. Contaminants go before Stage 02 because cycloess estimates its reference from each
# sample's own intensity distribution.

pacman::p_load(
  proteoDA, here, readxl, readr, dplyr, tidyr, tibble, stringr, purrr, openxlsx
)
source(here("04_Figures", "shared", "pca.R"))
source(here("01_Filtering", "a_script", "filter_config.R"))
cfg <- filter_cfg

data_dir <- here("01_Filtering", "c_data")
report_dir <- here("01_Filtering", "b_reports")
clear_dir <- function(d) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  unlink(setdiff(list.files(d, full.names = TRUE), file.path(d, ".gitkeep")), recursive = TRUE)
}
walk(c(data_dir, report_dir), clear_dir)

strip_iso <- function(x) sub("-[0-9]+$", "", x)

raw <- read_excel(here("00_input", "HRvLR_raw.xlsx"))
annot_cols <- c("uniprot_id", "protein", "gene", "description", "n_seq")
metadata <- as.data.frame(read_csv(here("00_input", "HRvLR_meta.csv"), show_col_types = FALSE))
rownames(metadata) <- metadata$Col_ID

annotation <- raw[, annot_cols]
intensity <- data.matrix(raw[, metadata$Col_ID])
stopifnot("Sample mismatch" = setequal(colnames(intensity), metadata$Col_ID))
n_raw <- nrow(annotation)

if (any(duplicated(annotation$uniprot_id))) { # keep the highest-mean row per accession
  keep_idx <- tibble(i = seq_len(n_raw), id = annotation$uniprot_id, m = rowMeans(intensity, na.rm = TRUE)) |>
    slice_max(m, n = 1, by = id, with_ties = FALSE) |>
    pull(i)
  annotation <- annotation[keep_idx, ]
  intensity <- intensity[keep_idx, ]
}

# Contamination is measured on the full matrix, before anything is removed. Deleting the
# evidence and then asserting the samples were clean proves nothing.
qc_index <- imap(qc_panels, \(genes, panel) {
  tibble(
    Col_ID = colnames(intensity),
    panel = panel,
    pct_signal = 100 * colSums(intensity[annotation$gene %in% genes, , drop = FALSE], na.rm = TRUE) /
      colSums(intensity, na.rm = TRUE)
  )
}) |>
  list_rbind() |>
  left_join(select(metadata, Col_ID, Subject_ID, Group, Timepoint), by = "Col_ID")

hpa <- read_tsv(here("00_input", cfg$hpa_file), show_col_types = FALSE) |>
  transmute(
    acc           = Uniprot,
    protein_class = `Protein class`,
    secretome     = `Secretome location`,
    blood_conc    = suppressWarnings(as.numeric(`Blood concentration - Conc. blood MS [pg/L]`)),
    ery           = suppressWarnings(as.numeric(`Single Cell Type RNA - Erythrocytes [nCPM]`)),
    myo           = suppressWarnings(as.numeric(`Single Cell Type RNA - Myonuclei [nCPM]`))
  ) |>
  filter(!is.na(acc), acc != "") |>
  separate_longer_delim(acc, delim = ", ") |>
  mutate(acc = strip_iso(acc)) |>
  distinct(acc, .keep_all = TRUE)

protein_calls <- annotation |>
  mutate(acc = strip_iso(uniprot_id)) |>
  left_join(hpa, by = "acc") |>
  left_join(select(contaminants, acc = uniprot_id, contam_class = class, contam_reason = reason),
    by = "acc"
  ) |>
  mutate(
    is_ery = !is.na(ery) & ery >= cfg$ery_cut,
    is_ig = str_detect(coalesce(protein_class, ""), "Immunoglobulin genes"),
    is_plasma = !is.na(secretome) & secretome == "Secreted to blood",
    rescued = !is.na(myo) & myo >= cfg$myo_cut & (is.na(blood_conc) | blood_conc < cfg$blood_max),
    verdict = case_when(
      !is.na(contam_class) ~ paste0("remove: ", contam_class),
      (is_ery | is_ig | is_plasma) & rescued ~ "keep: rescued (muscle-expressed)",
      is_ery ~ "remove: erythrocyte",
      is_ig ~ "remove: immunoglobulin",
      is_plasma ~ "remove: plasma",
      TRUE ~ "keep"
    ),
    reason = coalesce(contam_reason, verdict)
  ) |>
  select(uniprot_id, gene, description, verdict, reason, secretome, blood_conc, ery, myo)

keep <- !str_starts(protein_calls$verdict, "remove")

int_df <- as.data.frame(intensity[keep, ])
annot_df <- as.data.frame(annotation[keep, ])
rownames(int_df) <- rownames(annot_df) <- annot_df$uniprot_id
dal <- zero_to_missing(DAList(data = int_df, annotation = annot_df, metadata = metadata))

# Outlier consensus: 4 methods, >=outlier_k must agree.
lg <- log2(dal$data + 1)
pct_missing <- colMeans(is.na(dal$data)) * 100
samp_med <- apply(lg, 2, median, na.rm = TRUE)
med_cor <- apply(cor(lg, use = "pairwise.complete.obs"), 2, \(x) median(x[x < 1], na.rm = TRUE))
pc3 <- run_pca(dal$data[rowSums(is.na(dal$data)) == 0, ], dal$metadata, log_transform = TRUE)$pca$x[, 1:3]

hampel <- function(x, k = cfg$mad_k) abs(x - median(x)) > k * mad(x)
tukey <- function(x) x > quantile(x, 0.75) + 1.5 * IQR(x)

outlier_diag <- dal$metadata |>
  select(Col_ID, Subject_ID, Group, Timepoint) |>
  mutate(
    pct_missing = pct_missing[Col_ID],
    delta_missing = ave(pct_missing, str_remove(Col_ID, "_T[123]$"),
      FUN = \(x) sapply(x, \(v) max(abs(v - x)))
    ),
    miss_flag = tukey(pct_missing) | tukey(delta_missing),
    pca_flag = mahalanobis(pc3, colMeans(pc3), cov(pc3)) > qchisq(1 - cfg$mahal_p, df = 3),
    mad_flag = hampel(samp_med[Col_ID]),
    cor_flag = med_cor[Col_ID] < median(med_cor) - cfg$mad_k * mad(med_cor),
    n_flags = miss_flag + pca_flag + mad_flag + cor_flag,
    consensus_outlier = n_flags >= cfg$outlier_k
  )
outlier_ids <- outlier_diag$Col_ID[outlier_diag$consensus_outlier]
if (length(outlier_ids)) dal <- filter_samples(dal, !(Col_ID %in% outlier_ids))

# Missingness runs last, on the samples that survive. Running it first counted detections from
# samples that were about to be discarded.
n_pre_miss <- nrow(dal$data)
dal <- filter_proteins_by_group(dal,
  min_reps = cfg$miss_min_reps,
  min_groups = cfg$miss_min_groups, grouping_column = "Group_Time"
)

removed <- count(filter(protein_calls, str_starts(verdict, "remove")), verdict, name = "n_removed")
filter_log <- tibble(
  step = c("Raw input", removed$verdict, "Outlier samples", "Missingness"),
  n_removed = c(NA_integer_, removed$n_removed, 0L, n_pre_miss - nrow(dal$data))
) |>
  mutate(
    n_after = n_raw - cumsum(coalesce(n_removed, 0L)),
    pct_of_raw = round(n_after / n_raw * 100, 1)
  )
print(as.data.frame(filter_log))

cat("\nContamination (% of sample signal):\n")
qc_index |>
  summarise(median = median(pct_signal), max = max(pct_signal), .by = panel) |>
  arrange(desc(median)) |>
  as.data.frame() |>
  print(digits = 3)

saveRDS(dal, file.path(data_dir, "DAList_filtered.rds"))
write_csv(qc_index, file.path(data_dir, "contamination_index.csv"))
write_csv(protein_calls, file.path(data_dir, "protein_calls.csv"))

write.xlsx(
  list(
    filter_log = filter_log,
    removed = filter(protein_calls, str_starts(verdict, "remove")) |> arrange(verdict, gene),
    rescued = filter(protein_calls, str_detect(verdict, "rescued")) |> arrange(desc(myo)),
    contamination_index = qc_index,
    outlier_diagnostics = outlier_diag
  ),
  file.path(data_dir, "filtering_report.xlsx"),
  overwrite = TRUE
)

saveRDS(
  list(
    cfg = cfg, n_raw = n_raw, filter_log = filter_log, protein_calls = protein_calls,
    qc_index = qc_index, outlier_diag = outlier_diag, outlier_ids = outlier_ids
  ),
  file.path(data_dir, "filtering_intermediates.rds")
)

cat(sprintf("\nDone: %d proteins x %d samples\n", nrow(dal$data), ncol(dal$data)))

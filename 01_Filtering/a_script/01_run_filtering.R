#!/usr/bin/env Rscript
# HRvLR Stage 01: UniProt dedup -> HPA presence filter -> myonuclei-rescue blood
# removal -> group-wise missingness -> 4-method outlier consensus. Writes the
# filtered, un-normalized DAList handed to Stage 02.
#
# Blood logic follows CvH: a protein is removed iff blood-derived (secreted-to-blood |
# immunoglobulin | erythrocyte-high) AND NOT myonuclei-expressed (rescued as real muscle).
# Source: Human Protein Atlas (Uhlen 2015 Science; 2019 blood atlas).

pacman::p_load(
  proteoDA, here, readxl, readr, dplyr, tidyr, stringr, openxlsx,
  ggplot2, forcats, patchwork
)
source(here("shared", "pca.R"))

data_dir <- here("01_Filtering", "c_data")
report_dir <- here("01_Filtering", "b_reports")
clear_dir <- function(d) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  unlink(setdiff(list.files(d, full.names = TRUE), file.path(d, ".gitkeep")), recursive = TRUE)
}
clear_dir(data_dir)
clear_dir(report_dir)

cfg <- list(
  miss_min_reps   = 5, # min detected samples per group (of 8)
  miss_min_groups = 1, # min Group_Time levels passing the threshold
  outlier_k       = 3, # methods that must agree for consensus (>=3/4)
  mahal_p         = 0.01, # PCA Mahalanobis chi-sq cutoff
  mad_k           = 3, # MAD multiplier for median intensity / correlation
  ery_cut         = 5000, # erythrocyte nCPM at/above = hemoglobin-class red-cell protein
  myo_cut         = 50, # myonuclei nCPM at/above = candidate for muscle rescue
  blood_max       = 1e9 # rescue only if blood conc below this; above = true plasma protein
)

# Load matrix + metadata

raw <- read_excel(here("00_input", "HRvLR_raw.xlsx"))
annot_cols <- c("uniprot_id", "protein", "gene", "description", "n_seq")
annotation <- raw[, annot_cols]
intensity <- raw[, setdiff(names(raw), annot_cols)]

metadata <- as.data.frame(read_csv(here("00_input", "HRvLR_meta.csv"), show_col_types = FALSE))
rownames(metadata) <- metadata$Col_ID
stopifnot("Sample mismatch" = setequal(colnames(intensity), metadata$Col_ID))
intensity <- intensity[, metadata$Col_ID]
n_raw <- nrow(annotation)
cat(sprintf("Raw: %d proteins x %d samples\n", n_raw, ncol(intensity)))

if (any(duplicated(annotation$uniprot_id))) { # keep highest-mean row per accession
  rm_mean <- rowMeans(data.matrix(intensity), na.rm = TRUE)
  keep_idx <- tibble(i = seq_along(rm_mean), id = annotation$uniprot_id, m = rm_mean) |>
    group_by(id) |>
    slice_max(m, n = 1, with_ties = FALSE) |>
    pull(i)
  annotation <- annotation[keep_idx, ]
  intensity <- intensity[keep_idx, ]
  cat(sprintf("Deduplicated: %d proteins\n", nrow(annotation)))
}

# HPA presence + myonuclei-rescue blood removal

bl <- read_tsv(here("00_input", "HPA_annotations.tsv"), show_col_types = FALSE) |>
  transmute(
    gene          = Gene,
    uniprot       = Uniprot,
    protein_class = `Protein class`,
    secretome     = `Secretome location`,
    blood_conc    = suppressWarnings(as.numeric(`Blood concentration - Conc. blood MS [pg/L]`)),
    ery           = suppressWarnings(as.numeric(`Single Cell Type RNA - Erythrocytes [nCPM]`)),
    myo           = suppressWarnings(as.numeric(`Single Cell Type RNA - Myonuclei [nCPM]`))
  ) |>
  mutate(
    secreted_blood = secretome == "Secreted to blood" & !is.na(secretome),
    is_ig = str_detect(coalesce(protein_class, ""), "Immunoglobulin genes"),
    is_erythrocyte = !is.na(ery) & ery >= cfg$ery_cut,
    muscle_keep = !is.na(myo) & myo >= cfg$myo_cut &
      (is.na(blood_conc) | blood_conc < cfg$blood_max),
    reason = case_when(
      muscle_keep & (secreted_blood | is_ig | is_erythrocyte) ~ "keep: muscle-expressed (rescued)",
      is_erythrocyte ~ "remove: erythrocyte (hemoglobin)",
      is_ig ~ "remove: immunoglobulin",
      secreted_blood ~ "remove: secreted-to-blood (plasma)",
      TRUE ~ "keep: not blood-associated"
    ),
    verdict = if_else(str_starts(reason, "remove"), "remove", "keep")
  )
strip_iso <- function(x) sub("-\\d+$", "", x)
hpa_acc <- bl |>
  filter(!is.na(uniprot), uniprot != "") |>
  tidyr::separate_longer_delim(uniprot, delim = ", ") |>
  mutate(acc = strip_iso(uniprot))
hpa_present <- unique(hpa_acc$acc)
remove_acc <- unique(hpa_acc$acc[hpa_acc$verdict == "remove"])

flog <- tibble(step = "Raw input", n_after = nrow(annotation), n_removed = NA_integer_)
acc <- strip_iso(annotation$uniprot_id)
keep <- acc %in% hpa_present
dropped_ids <- annotation$uniprot_id[!keep]
flog <- bind_rows(flog, tibble(step = "HPA presence", n_after = sum(keep), n_removed = sum(!keep)))
annotation <- annotation[keep, ]
intensity <- intensity[keep, ]
readr::write_csv(
  tibble(uniprot_id = dropped_ids),
  here("01_Filtering", "c_data", "hpa_absent_dropped.csv")
)
acc <- strip_iso(annotation$uniprot_id)
keep <- !(acc %in% remove_acc)
flog <- bind_rows(flog, tibble(step = "Blood contaminant removal", n_after = sum(keep), n_removed = sum(!keep)))
annotation <- annotation[keep, ]
intensity <- intensity[keep, ]

# Build DAList + missingness filter

int_mat <- as.data.frame(data.matrix(intensity))
rownames(int_mat) <- annotation$uniprot_id
annot_df <- as.data.frame(annotation)
rownames(annot_df) <- annotation$uniprot_id
dal <- zero_to_missing(DAList(data = int_mat, annotation = annot_df, metadata = metadata))

n0 <- nrow(dal$data)
dal <- filter_proteins_by_group(dal,
  min_reps = cfg$miss_min_reps,
  min_groups = cfg$miss_min_groups, grouping_column = "Group_Time"
)
flog <- bind_rows(flog, tibble(
  step = sprintf("Missingness (>=%d of 8 in >=%d group)", cfg$miss_min_reps, cfg$miss_min_groups),
  n_after = nrow(dal$data), n_removed = n0 - nrow(dal$data)
))
flog <- flog |> mutate(pct_of_raw = round(n_after / n_raw * 100, 1))
print(as.data.frame(flog))

filtered_proteins <- annot_df |>
  filter(!uniprot_id %in% rownames(dal$data)) |>
  select(uniprot_id, gene, description)

# Outlier consensus (4-method, >=3/4)

# Method 1: sample missingness, with a within-subject delta arm
pct_missing <- colMeans(is.na(dal$data)) * 100
miss_info <- dal$metadata |>
  select(Col_ID, Subject_ID, Group, Timepoint) |>
  mutate(pct_missing = pct_missing[Col_ID], prefix = str_remove(Col_ID, "_T[123]$")) |>
  group_by(prefix) |>
  mutate(delta_missing = sapply(pct_missing, \(x) max(abs(x - pct_missing)))) |>
  ungroup()
miss_thresh <- quantile(pct_missing, 0.75) + 1.5 * IQR(pct_missing)
delta_thresh <- quantile(miss_info$delta_missing, 0.75) + 1.5 * IQR(miss_info$delta_missing)
miss_info$miss_flag <- miss_info$pct_missing > miss_thresh | miss_info$delta_missing > delta_thresh

# Method 2: PCA Mahalanobis distance on PC1-3
complete_mat <- dal$data[rowSums(is.na(dal$data)) == 0, ]
cat(sprintf("PCA on %d complete proteins (of %d)\n", nrow(complete_mat), nrow(dal$data)))
pca_pre <- run_pca(complete_mat, dal$metadata, log_transform = TRUE)
pc3 <- pca_pre$pca$x[, 1:3]
mahal <- mahalanobis(pc3, colMeans(pc3), cov(pc3))
pca_flags <- tibble(
  Col_ID = colnames(dal$data), mahal_dist = mahal,
  pca_flag = mahal > qchisq(1 - cfg$mahal_p, df = 3)
)

# Method 3: MAD-based median intensity
samp_med <- apply(log2(dal$data + 1), 2, median, na.rm = TRUE)
global_med <- median(samp_med)
mad_val <- mad(samp_med)
mad_flags <- tibble(
  Col_ID = names(samp_med), sample_median = samp_med,
  mad_flag = abs(samp_med - global_med) > cfg$mad_k * mad_val
)

# Method 4: inter-sample correlation
cor_mat <- cor(log2(dal$data + 1), use = "pairwise.complete.obs")
med_cor <- apply(cor_mat, 2, function(x) median(x[x < 1], na.rm = TRUE))
cor_med_all <- median(med_cor)
cor_mad_all <- mad(med_cor)
cor_flags <- tibble(
  Col_ID = names(med_cor), median_cor = med_cor,
  cor_flag = med_cor < cor_med_all - cfg$mad_k * cor_mad_all
)

outlier_diag <- miss_info |>
  left_join(pca_flags, by = "Col_ID") |>
  left_join(mad_flags, by = "Col_ID") |>
  left_join(cor_flags, by = "Col_ID") |>
  mutate(
    n_flags = miss_flag + pca_flag + mad_flag + cor_flag,
    consensus_outlier = n_flags >= cfg$outlier_k
  )
outlier_ids <- outlier_diag$Col_ID[outlier_diag$consensus_outlier]
n_outliers <- length(outlier_ids)
cat(sprintf(
  "Outliers (>=%d/4): %s\n", cfg$outlier_k,
  if (n_outliers) paste(outlier_ids, collapse = ", ") else "none"
))

data_pre_outlier <- dal$data
meta_pre_outlier <- dal$metadata
if (n_outliers > 0) dal <- filter_samples(dal, !(Col_ID %in% outlier_ids))
cat(sprintf("%d samples remain\n", ncol(dal$data)))

# Export

saveRDS(dal, file.path(data_dir, "DAList_filtered.rds"))

ann_lookup <- raw |>
  as.data.frame() |>
  distinct(gene, .keep_all = TRUE) |>
  select(gene, protein, description)
removed <- bl |>
  filter(verdict == "remove", gene %in% raw$gene) |>
  left_join(ann_lookup, by = "gene") |>
  transmute(gene, uniprot, protein, description, reason, secretome,
    blood_conc,
    erythrocyte_nCPM = ery, myonuclei_nCPM = myo
  ) |>
  arrange(reason, desc(blood_conc))
rescued <- bl |>
  filter(str_detect(reason, "rescued"), gene %in% raw$gene) |>
  left_join(ann_lookup, by = "gene") |>
  transmute(gene, uniprot, protein, description, reason, secretome,
    blood_conc,
    erythrocyte_nCPM = ery, myonuclei_nCPM = myo
  ) |>
  arrange(desc(myonuclei_nCPM))

write.xlsx(
  list(
    filter_log = flog,
    contaminants_removed = removed,
    contaminants_rescued = rescued,
    outlier_diagnostics = outlier_diag,
    blood_classification = bl |> select(
      gene, uniprot, protein_class, secretome,
      blood_conc, ery, myo, reason, verdict
    )
  ),
  file.path(data_dir, "filtering_report.xlsx"),
  overwrite = TRUE
)

pal <- c(
  "remove: secreted-to-blood (plasma)" = "#D6604D", "remove: immunoglobulin" = "#F4A582",
  "remove: erythrocyte (hemoglobin)" = "#9970AB"
)
bars <-
  (count(removed, reason) |> ggplot(aes(reorder(reason, n), n, fill = reason)) +
    geom_col(show.legend = FALSE) +
    geom_text(aes(label = n), hjust = -0.2, size = 3.5) +
    scale_fill_manual(values = pal) +
    coord_flip(clip = "off") +
    labs(title = "A  Filtered out, by reason", x = NULL, y = "proteins")) /
    (removed |> filter(!is.na(blood_conc)) |> slice_max(blood_conc, n = 25) |>
      ggplot(aes(fct_reorder(gene, blood_conc), log10(blood_conc + 1), fill = reason)) +
      geom_col(show.legend = FALSE) +
      scale_fill_manual(values = pal) +
      coord_flip() +
      labs(
        title = "B  Top filtered-out plasma proteins by blood concentration",
        x = NULL, y = "log10 blood conc [pg/L]"
      )) +
    plot_layout(heights = c(1, 3)) & theme_minimal(base_size = 11)
ggsave(file.path(report_dir, "blood_contaminants.pdf"), bars, width = 8, height = 9)

filter_bar_data <- flog |>
  filter(!is.na(n_removed)) |>
  mutate(step = factor(step, levels = step)) |>
  pivot_longer(c(n_after, n_removed), names_to = "status", values_to = "n") |>
  mutate(status = recode(status, n_after = "Retained", n_removed = "Removed"))

miss_bar_data <- meta_pre_outlier |>
  select(Col_ID, Group_Time) |>
  mutate(
    detected = colSums(!is.na(data_pre_outlier[, Col_ID])),
    missing = nrow(data_pre_outlier) - detected,
    is_outlier = Col_ID %in% outlier_ids
  ) |>
  pivot_longer(c(detected, missing), names_to = "status", values_to = "n") |>
  mutate(status = str_to_title(status))

saveRDS(
  list(
    cfg = cfg, filter_log = flog, filter_bar_data = filter_bar_data,
    miss_bar_data = miss_bar_data, n_raw = n_raw, n_outliers = n_outliers,
    outlier_diag = outlier_diag, outlier_ids = outlier_ids,
    miss_thresh = miss_thresh, delta_thresh = delta_thresh, pca_pre = pca_pre,
    global_med = global_med, mad_val = mad_val, cor_med_all = cor_med_all,
    cor_mad_all = cor_mad_all, filtered_proteins = filtered_proteins,
    data_pre_outlier = data_pre_outlier, meta_pre_outlier = meta_pre_outlier
  ),
  file.path(data_dir, "filtering_intermediates.rds")
)

if (file.exists("Rplots.pdf")) file.remove("Rplots.pdf")
cat(sprintf("Done: %d proteins x %d samples -> %s/\n", nrow(dal$data), ncol(dal$data), data_dir))

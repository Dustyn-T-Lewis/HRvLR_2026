# Protein report to a filtered, unnormalised DAList. Contaminants go by identity, outlier samples
# by a four-method consensus, then proteins by detection on the samples that remain. No threshold
# reads a contrast, phenotype or group label.

suppressPackageStartupMessages({
  library(here)
  library(proteoDA)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(purrr)
  library(readr)
  library(readxl)
  library(ggplot2)
  library(forcats)
  library(patchwork)
  library(writexl)
})

out <- here("01_Preprocess", "01_Filtering", "c_data")
figures <- here("01_Preprocess", "01_Filtering", "b_reports")
walk(c(out, figures), dir.create, recursive = TRUE, showWarnings = FALSE)
inputs <- c(
  report = "00_Input/HRvLR_raw.xlsx",
  samples = "00_Input/metadata.csv",
  blood = "00_Input/blood_contaminants.csv",
  hpa = "00_Input/HPA_annotations_full.tsv",
  rbc = "00_Input/RBC_proteome_reference.tsv"
)
paths <- map_chr(inputs, here)
stopifnot(file.exists(paths))

strip_iso <- function(x) sub("-[0-9]+$", "", x)

cfg <- list(
  miss_min_reps = 5, # detected samples a group cell needs
  miss_min_groups = 1, # cells that must clear miss_min_reps
  cor_min_obs = 15, # samples a protein needs before blood_cor is computed
  outlier_k = 3, # outlier methods that must agree, of four
  mahal_p = 0.01, # PCA Mahalanobis chi-square tail
  mad_k = 3, # Hampel constant
  ery_cut = 5000, # erythrocyte nCPM at or above: red-cell protein
  # myonuclei nCPM at or above: candidate for muscle rescue. At 50 the rescue missed GPI (43.1),
  # ANXA2 (40.2) and PPIA (21.1); blood_max, not this floor, separates plasma from muscle.
  myo_cut = 20,
  blood_max = 1e9, # rescue only below this plasma concentration, pg/L
  blood_anchor = c("HBB", "HBA1", "HBD", "HBG1", "HBG2")
)

# Matched by accession only. cRAP collides with muscle on symbol (ALDOA against ALDOA_RABIT) and on
# description (a trypsin or albumin pattern takes parvalbumin). Only cRAP's dust and contact
# section is used: its UPS1 section holds myoglobin and creatine kinase M.
contaminants <- tribble(
  ~uniprot_id, ~gene, ~class, ~reason,
  "P04264", "KRT1", "keratin", "cRAP dust/contact; cornified epidermal",
  "P35908", "KRT2", "keratin", "cRAP dust/contact; epidermal",
  "P13645", "KRT10", "keratin", "cRAP dust/contact; type-I partner of KRT1",
  "P35527", "KRT9", "keratin", "cRAP dust/contact; palmoplantar",
  "P08779", "KRT16", "keratin", "cRAP dust/contact",
  "P04259", "KRT6B", "keratin", "cRAP dust/contact",
  "P02533", "KRT14", "keratin", "cRAP dust/contact; basal epidermal",
  "P13647", "KRT5", "keratin", "cRAP dust/contact; basal epidermal",
  "P02538", "KRT6A", "keratin", "cRAP dust/contact",
  "Q04695", "KRT17", "keratin", "cRAP dust/contact",
  "Q7Z794", "KRT77", "keratin", "cRAP dust/contact",
  "Q8N1N4", "KRT78", "keratin", "cRAP dust/contact",
  "P05787", "KRT8", "keratin", "simple-epithelial; Spearman +0.78 with the KRT1/2/10 load",
  "P15924", "DSP", "keratin", "desmoplakin; skin desmosome",
  "Q02413", "DSG1", "keratin", "desmoglein-1; skin desmosome",
  "P81605", "DCD", "keratin", "dermcidin; sweat",
  "P31151", "S100A7", "keratin", "psoriasin; skin",
  "P01040", "CSTA", "keratin", "cystatin-A; cornified envelope",
  "P68871", "HBB", "globin", "red-cell carryover",
  "P69905", "HBA1", "globin", "red-cell carryover",
  "P02042", "HBD", "globin", "red-cell carryover",
  "P69891", "HBG1", "globin", "red-cell carryover",
  "P69892", "HBG2", "globin", "orphan gamma-globin; shares tryptic peptides with HBB/HBD",
  "P09105", "HBQ1", "globin", "red-cell carryover",
  "Q9NZD4", "AHSP", "globin", "erythroid-specific haemoglobin chaperone",
  "P01834", "IGKC", "immunoglobulin", "absent from HPA; constant region",
  "P01860", "IGHG3", "immunoglobulin", "absent from HPA; constant region",
  "P01880", "IGHD", "immunoglobulin", "absent from HPA; constant region",
  "P02746", "C1QB", "complement", "absent from HPA",
  "P00736", "C1R", "complement", "absent from HPA"
) |>
  bind_rows(
    read_csv(paths[["blood"]], col_types = cols(.default = "c")) |>
      select(uniprot_id, gene, class, reason)
  )

# Share of each sample's summed intensity, measured before removal. Descriptive only.
qc_panels <- list(
  keratin = c(
    "KRT1", "KRT2", "KRT10", "KRT8", "KRT9", "KRT16", "KRT14", "KRT5", "KRT6A",
    "KRT6B", "DCD"
  ),
  blood = c(
    "HBB", "HBA1", "HBD", "HBG1", "HBG2", "HBQ1", "AHSP", "ALB", "CA1", "SLC4A1",
    "SPTA1"
  ),
  leukocyte = c("MPO", "CTSG", "ELANE", "PRTN3", "AZU1", "S100A8", "S100A9", "LYZ", "LTF"),
  adipose = c("FABP4", "PLIN1", "PLIN4", "ADIPOQ", "ADIRF"),
  myofibre = c(
    "MB", "CKM", "ACTA1", "MYH1", "MYH2", "MYH7", "ALDOA", "CASQ1", "PYGM",
    "ATP2A1", "DES"
  )
)

raw <- read_excel(paths[["report"]])
metadata <- read_csv(paths[["samples"]], show_col_types = FALSE) |>
  select(sample_id, subject, arm, timepoint, group) |>
  as.data.frame()
rownames(metadata) <- metadata$sample_id

annotation_cols <- c("uniprot_id", "protein", "gene", "description", "n_seq")
# A run missing from either side would drop out silently or come back all-NA.
runs <- setdiff(names(raw), annotation_cols)
stopifnot("sample sheet and report disagree" = setequal(runs, metadata$sample_id))
annotation <- raw[, annotation_cols]
intensity <- data.matrix(raw[, metadata$sample_id])
n_raw <- nrow(annotation)

# One row per accession: the highest-mean row wins.
keep_row <- tibble(
  i = seq_len(n_raw), id = annotation$uniprot_id, m = rowMeans(intensity, na.rm = TRUE)
) |>
  slice_max(m, n = 1, by = id, with_ties = FALSE) |>
  pull(i)
annotation <- annotation[keep_row, ]
intensity <- intensity[keep_row, ]
n_dedup <- length(keep_row)

# The blood index is each sample's mean log2 haemoglobin, taken before the anchors are removed.
# blood_cor, each protein's Spearman correlation with it, is reported and gates nothing.
log_int <- log2(intensity)
log_int[!is.finite(log_int)] <- NA
blood_index <- colMeans(log_int[annotation$gene %in% cfg$blood_anchor, , drop = FALSE],
  na.rm = TRUE
)
testable <- rowSums(!is.na(log_int)) >= cfg$cor_min_obs
blood_cor <- suppressWarnings(as.vector(cor(
  t(log_int), blood_index,
  method = "spearman", use = "pairwise.complete.obs"
)))
blood_cor[!testable] <- NA

qc_index <- imap(qc_panels, \(genes, panel) {
  tibble(
    sample_id = colnames(intensity), panel = panel,
    pct_signal = 100 * colSums(intensity[annotation$gene %in% genes, , drop = FALSE],
      na.rm = TRUE
    ) / colSums(intensity, na.rm = TRUE)
  )
}) |>
  list_rbind() |>
  left_join(select(metadata, sample_id, subject, arm, timepoint), by = "sample_id")

# A protein goes when it is on the curated list, or when HPA marks it plasma, immunoglobulin or
# erythrocyte and the muscle rescue does not reach it. Absence from HPA never removes. Red-cell
# membership is reported as in_rbc and removes nothing.
hpa <- read_tsv(paths[["hpa"]], show_col_types = FALSE) |>
  transmute(
    acc = Uniprot, protein_class = `Protein class`, secretome = `Secretome location`,
    blood_conc = suppressWarnings(as.numeric(`Blood concentration - Conc. blood MS [pg/L]`)),
    ery = suppressWarnings(as.numeric(`Single Cell Type RNA - Erythrocytes [nCPM]`)),
    myo = suppressWarnings(as.numeric(`Single Cell Type RNA - Myonuclei [nCPM]`))
  ) |>
  filter(!is.na(acc), acc != "") |>
  separate_longer_delim(acc, delim = ", ") |>
  mutate(acc = strip_iso(acc)) |>
  distinct(acc, .keep_all = TRUE)
rbc <- read_tsv(paths[["rbc"]], show_col_types = FALSE)

protein_calls <- annotation |>
  mutate(acc = strip_iso(uniprot_id), blood_cor = blood_cor) |>
  left_join(hpa, by = "acc") |>
  left_join(
    select(contaminants, acc = uniprot_id, contam_class = class, contam_reason = reason),
    by = "acc"
  ) |>
  mutate(
    is_curated = !is.na(contam_class),
    is_ery = !is.na(ery) & ery >= cfg$ery_cut,
    is_ig = str_detect(coalesce(protein_class, ""), "Immunoglobulin genes"),
    is_plasma = !is.na(secretome) & secretome == "Secreted to blood",
    in_rbc = acc %in% na.omit(rbc$acc) | gene %in% na.omit(rbc$gene),
    rescued = !is.na(myo) & myo >= cfg$myo_cut &
      (is.na(blood_conc) | blood_conc < cfg$blood_max),
    is_blood = is_ery | is_ig | is_plasma,
    contaminant = is_curated | (is_blood & !rescued),
    verdict = case_when(
      is_curated ~ paste0("remove: ", contam_class),
      is_blood & rescued ~ "keep: rescued (muscle-expressed)",
      is_plasma ~ "remove: plasma",
      is_ig ~ "remove: immunoglobulin",
      is_ery ~ "remove: erythrocyte",
      TRUE ~ "keep"
    ),
    reason = coalesce(contam_reason, verdict)
  ) |>
  select(
    uniprot_id, gene, description, verdict, reason, contaminant, in_rbc,
    secretome, blood_conc, ery, myo, blood_cor
  )
stopifnot(protein_calls$contaminant == str_starts(protein_calls$verdict, "remove"))

# A plain subset: filter_proteins_by_annotation() calls if() on a length-2 class vector and errors
# on every real DAList.
keep <- !protein_calls$contaminant
int_df <- as.data.frame(intensity[keep, ])
annot_df <- as.data.frame(annotation[keep, ])
rownames(int_df) <- rownames(annot_df) <- annot_df$uniprot_id
dal <- zero_to_missing(DAList(data = int_df, annotation = annot_df, metadata = metadata))

# Four methods flag a sample: missingness (Tukey fence on its share and on its spread within
# subject), PCA Mahalanobis distance, median intensity (Hampel) and median inter-sample
# correlation. Three of four remove it.
lg <- log2(dal$data + 1)
complete <- dal$data[rowSums(is.na(dal$data)) == 0, ]
pcs <- prcomp(t(log2(complete + 1)), center = TRUE, scale. = TRUE)$x[, 1:3]
med_cor <- apply(cor(lg, use = "pairwise.complete.obs"), 2, \(x) median(x[x < 1], na.rm = TRUE))
hampel <- function(x) abs(x - median(x)) > cfg$mad_k * mad(x)
tukey <- function(x) x > quantile(x, 0.75) + 1.5 * IQR(x)

outlier_diag <- dal$metadata |>
  select(sample_id, subject, arm, timepoint) |>
  mutate(
    pct_missing = colMeans(is.na(dal$data))[sample_id] * 100,
    delta_missing = ave(pct_missing, subject, FUN = \(x) pmax(x - min(x), max(x) - x)),
    miss_flag = tukey(pct_missing) | tukey(delta_missing),
    pca_flag = mahalanobis(pcs, colMeans(pcs), cov(pcs)) > qchisq(1 - cfg$mahal_p, df = 3),
    mad_flag = hampel(apply(lg, 2, median, na.rm = TRUE)[sample_id]),
    cor_flag = med_cor[sample_id] < median(med_cor) - cfg$mad_k * mad(med_cor),
    n_flags = miss_flag + pca_flag + mad_flag + cor_flag,
    consensus_outlier = n_flags >= cfg$outlier_k
  )
stopifnot(identical(rownames(pcs), outlier_diag$sample_id))
outlier_ids <- outlier_diag$sample_id[outlier_diag$consensus_outlier]
dal <- filter_samples(dal, !(sample_id %in% outlier_ids))

# Detection runs last, so no discarded sample counts toward a protein's detections.
n_before <- nrow(dal$data)
dal <- filter_proteins_by_group(dal,
  min_reps = cfg$miss_min_reps, min_groups = cfg$miss_min_groups,
  grouping_column = "group"
)
removed <- count(filter(protein_calls, contaminant), verdict, name = "n_removed")
filter_log <- tibble(
  step = c("raw input", "duplicate accession", removed$verdict, "outlier samples", "detection"),
  n_removed = c(NA, n_raw - n_dedup, removed$n_removed, 0L, n_before - nrow(dal$data))
) |>
  mutate(
    n_after = n_raw - cumsum(coalesce(n_removed, 0L)),
    pct_of_raw = round(100 * n_after / n_raw, 1),
    n_samples = c(rep(nrow(metadata), nrow(removed) + 2), rep(ncol(dal$data), 2))
  )

set.seed(42) # the jittered points only
waterfall <- filter_log |>
  filter(!is.na(n_removed)) |>
  mutate(step = fct_inorder(str_remove(step, "^remove: ")), ymax = n_after + n_removed)
cascade <- ggplot(waterfall, aes(step)) +
  geom_rect(aes(
    xmin = as.numeric(step) - 0.35, xmax = as.numeric(step) + 0.35,
    ymin = n_after, ymax = ymax
  ), fill = "#D6604D") +
  geom_step(aes(y = n_after, group = 1), direction = "mid", colour = "grey35") +
  geom_text(aes(y = ymax, label = if_else(n_removed > 0, paste0("-", n_removed), "")),
    vjust = -0.6, size = 3, colour = "#B2182B"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(
    title = "Filtering cascade",
    subtitle = sprintf("%d proteins measured, %d analysed", n_raw, nrow(dal$data)),
    x = NULL, y = "proteins retained"
  ) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

contamination <- qc_index |>
  filter(panel != "myofibre") |>
  ggplot(aes(timepoint, pct_signal, fill = timepoint)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.75) +
  geom_jitter(aes(shape = arm), width = 0.15, size = 1.4) +
  facet_wrap(~panel, scales = "free_y", nrow = 1) +
  scale_fill_brewer(palette = "Blues", guide = "none") +
  scale_shape_manual(values = c(HR = 16, LR = 1), name = NULL) +
  labs(
    title = "Contamination per sample, before removal",
    subtitle = "share of summed intensity per panel", x = NULL, y = "% of sample signal"
  )

flags <- outlier_diag |>
  filter(n_flags > 0) |>
  pivot_longer(ends_with("_flag"), names_to = "method", values_to = "flagged") |>
  mutate(sample_id = fct_reorder(sample_id, n_flags))
consensus <- ggplot(flags, aes(method, sample_id, fill = flagged)) +
  geom_tile(colour = "white") +
  scale_fill_manual(values = c(`TRUE` = "#D6604D", `FALSE` = "grey92"), name = NULL) +
  labs(
    title = "Outlier flags, flagged samples only",
    subtitle = sprintf(
      "%d of %d samples flagged; %d met the %d-of-4 consensus",
      n_distinct(flags$sample_id), nrow(metadata), length(outlier_ids), cfg$outlier_k
    ),
    x = NULL, y = NULL
  )

filter_figure <- (cascade / contamination / consensus) +
  plot_layout(heights = c(1, 1.2, 0.9)) &
  theme_minimal(base_size = 10)

# Does blood content rise more at T3 in one arm? Mixed model on the analysed samples (arm by
# timepoint, subject intercept) on two scales; the permutation p shuffles the arm label across
# subjects 999 times. A confound on this scale does not cancel in the interaction contrasts.
blood_data <- dal$metadata |>
  select(sample_id, subject, arm, timepoint) |>
  mutate(
    arm = factor(arm, c("HR", "LR")),
    log2_index = unname(blood_index[sample_id]),
  ) |>
  left_join(
    filter(qc_index, panel == "blood") |> select(sample_id, pct_signal),
    by = "sample_id"
  )
arm_by_t3 <- function(data, response) {
  fit <- suppressMessages(lme4::lmer(
    reformulate(c("arm * timepoint", "(1 | subject)"), response),
    data = data
  ))
  coef(summary(fit))["armLR:timepointT3", c("Estimate", "t value")]
}
set.seed(42)
arms <- distinct(blood_data, subject, arm)
blood_by_arm <- map(set_names(c("log2_index", "pct_signal")), \(response) {
  observed <- arm_by_t3(blood_data, response)
  null_t <- replicate(999, {
    shuffled <- set_names(sample(arms$arm), arms$subject)
    arm_by_t3(mutate(blood_data, arm = shuffled[subject]), response)[["t value"]]
  })
  tibble(
    scale = response, lr_minus_hr_at_t3 = round(observed[["Estimate"]], 3),
    t = round(observed[["t value"]], 2),
    permutation_p = (1 + sum(abs(null_t) >= abs(observed[["t value"]]))) / 1000
  )
}) |>
  list_rbind()

# The old blood_cor cut of 0.45 rested on a null recorded at 0.43. Recomputed and split by how many
# samples saw a protein, the null differs by observation count, so no single cut serves all.
set.seed(42)
tested <- log_int[testable, ]
null_rho <- replicate(300, abs(as.vector(suppressWarnings(cor(
  t(tested), sample(blood_index),
  method = "spearman", use = "pairwise.complete.obs"
)))))
n_obs <- cut(rowSums(!is.na(tested)), c(0, 25, 40, Inf), c("15-25", "26-40", "41-48"))
null_999 <- \(rows) round(quantile(null_rho[rows, ], 0.999, na.rm = TRUE, names = FALSE), 3)
blood_null <- bind_rows(
  map(levels(n_obs), \(bin) {
    tibble(observations = bin, proteins = sum(n_obs == bin), null_999 = null_999(n_obs == bin))
  }),
  tibble(observations = "all", proteins = length(n_obs), null_999 = null_999(TRUE))
)

saveRDS(dal, file.path(out, "DAList_filtered.rds"), compress = "xz")
pdf(file.path(figures, "01_filter_figures.pdf"), width = 11, height = 8.5)
print(filter_figure)
invisible(dev.off())
sheets <- list(
  filter_log = filter_log,
  contaminants_removed = arrange(filter(protein_calls, contaminant), verdict, gene),
  rescued = arrange(filter(protein_calls, str_detect(verdict, "rescued")), desc(myo)),
  protein_calls = protein_calls,
  contamination_index = qc_index,
  blood_index = tibble(sample_id = names(blood_index), blood_index = unname(blood_index)),
  outlier_diagnostics = outlier_diag,
  blood_by_arm = blood_by_arm,
  blood_cor_null = blood_null,
  input_manifest = tibble(
    input = names(inputs), path = unname(inputs), md5 = unname(tools::md5sum(paths))
  ),
  package_versions = sessioninfo::package_info("loaded", dependencies = FALSE) |>
    as_tibble() |>
    select(package, version = loadedversion, source)
)
read_me <- tibble(sheet = names(sheets), holds = c(
  "Proteins left after each step.",
  "Every protein removed as a contaminant, with its verdict and reason.",
  "Blood-tagged proteins kept by the muscle rescue.",
  "One verdict per protein, with the annotation behind it.",
  "Share of each sample's signal per contamination panel, before removal.",
  "Mean log2 haemoglobin intensity per sample.",
  "The four outlier flags per sample.",
  "Arm by T3 blood term on two scales, with a permutation p.",
  "99.9th percentile of the blood_cor permutation null, by observation count.",
  "Files read, with md5.",
  "Packages loaded at run time."
))
write_xlsx(c(list(read_me = read_me), sheets), file.path(out, "01_filter.xlsx"))

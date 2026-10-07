# limpa on the normalised protein matrix. A missing value says the protein was probably faint, so
# limpa fits one detection probability curve (DPC) for the run and gives every protein a value in
# every sample with a standard error: small where observed, wide where missing. Our input is
# already summarised to proteins, so each row is its own feature (dpcQuantByRow, the route the
# limpa vignette gives for protein-level input).
#
# Writes the limpa object for the protein fit, a complete matrix for fry, and a blood-corrected
# matrix for the methods that cannot carry a covariate.

suppressPackageStartupMessages({
  library(here)
  library(limma)
  library(limpa)
  library(dplyr)
  library(tibble)
})
source(here("R", "helpers.R"))

stage <- here("01_Processing", "03_Quantification")
out <- file.path(stage, "c_data")
figures <- file.path(stage, "b_reports")
for (path in c(out, figures)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
inputs <- c(normalized = "01_Processing/02_Normalisation/c_data/DAList_normalized.rds")
paths <- vapply(inputs, here, character(1))
stopifnot(file.exists(paths))
dal <- readRDS(paths[["normalized"]])
values <- as.matrix(dal$data)
stopifnot(identical(dal$metadata$sample_id, colnames(values)))

group_levels <- c("HR_T1", "HR_T2", "HR_T3", "LR_T1", "LR_T2", "LR_T3")
meta <- as_tibble(dal$metadata) |>
  mutate(
    group = factor(group, group_levels),
    blood = blood_rel - mean(blood_rel)
  )

curve <- dpc(values)
proteins <- dpcQuantByRow(values, dpc = curve, verbose = FALSE)
observed <- proteins$other$n.observations > 0
stopifnot(!anyNA(proteins$E), identical(observed, !is.na(values)))

# The blood slope is fitted beside time point only, so no group label shapes the corrected matrix
# that every analysis reads.
corrected <- blood_corrected(proteins$E, meta$blood, model.matrix(~ 0 + timepoint, meta))
as_dalist <- function(data) {
  list(data = data, annotation = dal$annotation, metadata = meta, observed = observed)
}
saveRDS(
  list(proteins = proteins, metadata = meta, annotation = dal$annotation),
  file.path(out, "proteins.rds"),
  compress = "xz"
)
saveRDS(as_dalist(proteins$E), file.path(out, "DAList_limpa.rds"), compress = "xz")
saveRDS(as_dalist(corrected), file.path(out, "DAList_corrected.rds"), compress = "xz")

pdf(file.path(figures, "quantification_report.pdf"), width = 11, height = 8.5)
plotDPC(curve)
title(sub = sprintf(
  "Detection probability curve: intercept %.2f, slope %.2f. %.1f%% of values missing.",
  curve$dpc[1], curve$dpc[2], 100 * mean(!observed)
))
invisible(dev.off())

standard_error <- proteins$other$standard.error
write_workbook(
  list(
    curve = tibble(
      term = names(curve$dpc), estimate = unname(curve$dpc),
      pct_missing = 100 * mean(!observed)
    ),
    proteins = tibble(
      uniprot_id = rownames(values), gene = dal$annotation$gene,
      n_observed = rowSums(observed), prop_observed = rowMeans(observed),
      mean_observed = rowMeans(values, na.rm = TRUE), mean_complete = rowMeans(proteins$E),
      median_se_observed = apply(
        ifelse(observed, standard_error, NA), 1, median,
        na.rm = TRUE
      ),
      median_se_missing = apply(ifelse(observed, NA, standard_error), 1, median, na.rm = TRUE)
    ),
    blood = select(meta, sample_id, subject, group, blood_rel, blood),
    # Pearson correlation between samples after centring each protein on its mean, so it reflects
    # how samples differ rather than the abundance order every sample shares.
    sample_correlation = as_tibble(cor(corrected - rowMeans(corrected)), rownames = "sample_id")
  ),
  c(
    "DPC intercept and slope on the logit scale, and the share of values missing.",
    paste(
      "Per protein: values observed, mean log2 abundance where observed and over the complete",
      "matrix, and the median standard error where observed and where missing."
    ),
    "Per sample: blood index (HBB and HBA1 against the sample median) and its centred value.",
    "Pearson correlation between samples, each protein centred, on the blood-corrected matrix."
  ),
  file.path(out, "quantification.xlsx"),
  inputs, paths
)

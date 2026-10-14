# Cyclic loess on the filtered DAList. Missing values stay missing: 03_Quantification models them
# with limpa, and blood enters every model as a covariate.

suppressPackageStartupMessages({
  library(here)
  library(proteoDA)
  library(dplyr)
  library(tibble)
  library(limma)
  library(ggplot2)
})
source(here("00_Setup", "Shared_functions", "analysis_reporting_helpers.R"))

inputs <- c(filtered = "01_Processing/01_Protein_Filtering/c_data/DAList_filtered.rds")
paths <- vapply(inputs, here, character(1))
stopifnot(file.exists(paths))
dal <- readRDS(paths[["filtered"]])
out <- here("01_Processing", "02_Normalisation", "c_data")
reports <- here("01_Processing", "02_Normalisation", "b_reports")
for (path in c(out, reports)) dir.create(path, recursive = TRUE, showWarnings = FALSE)

# proteoDA draws every method it offers, and QC reports before and after cyclic loess, into
# b_reports/. Cyclic loess is limma's fast method with an adaptive span (chooseLowessSpan), so the
# span moves when filtering changes the protein count. A null device catches the plot they draw to
# the default device, which would otherwise land in Rplots.pdf.
pdf(NULL)
write_norm_report(dal,
  grouping_column = "group", output_dir = reports,
  filename = "norm_comparison.pdf", overwrite = TRUE
)
write_qc_report(dal,
  color_column = "group", output_dir = reports,
  filename = "qc_pre.pdf", overwrite = TRUE
)
# normalize_data() logs the raw intensities, so the "before" values are logged here to match.
before <- log2(as.matrix(dal$data))
dal <- normalize_data(dal, norm_method = "cycloess")
write_qc_report(dal,
  color_column = "group", output_dir = reports,
  filename = "qc_post.pdf", overwrite = TRUE
)
invisible(dev.off())

stopifnot(identical(dal$metadata$sample_id, colnames(dal$data)))
# eta squared: share of each protein's variance the six group cells explain.
cell_fit <- lmFit(dal$data, model.matrix(~group, dal$metadata))
total_ss <- apply(dal$data, 1, \(x) sum((x - mean(x, na.rm = TRUE))^2, na.rm = TRUE))
eta2 <- 1 - cell_fit$sigma^2 * cell_fit$df.residual / total_ss
eta2[rowSums(!is.na(dal$data)) < 4] <- NA
filled <- apply(dal$data, 2, \(x) replace(x, is.na(x), median(x, na.rm = TRUE)))
pca <- prcomp(t(filled), center = TRUE, scale. = TRUE)
components <- tibble(
  component = paste0("PC", 1:4),
  variance_explained = round(summary(pca)$importance[2, 1:4], 3)
)

pca_scores <- as_tibble(pca$x[, 1:2], rownames = "sample_id") |>
  left_join(dal$metadata, by = "sample_id")
pca_figure <- ggplot(pca_scores, aes(PC1, PC2, colour = timepoint, shape = arm)) +
  geom_point(size = 2.2) +
  scale_colour_manual(values = c(T1 = "#E69F00", T2 = "#0072B2", T3 = "#009E73")) +
  labs(
    title = "Samples after normalisation",
    subtitle = sprintf("PCA of %d proteins, gaps median-filled for display", nrow(filled)),
    x = sprintf("PC1 (%.1f%%)", 100 * components$variance_explained[1]),
    y = sprintf("PC2 (%.1f%%)", 100 * components$variance_explained[2])
  ) +
  theme_minimal(base_size = 10)

saveRDS(dal, file.path(out, "DAList_normalized.rds"), compress = "xz")
pdf(file.path(reports, "normalisation_report.pdf"), width = 11, height = 8.5)
print(pca_figure)
invisible(dev.off())
# Each sample's log2 intensity quantiles before and after cyclic loess, for the QC figure.
intensity_quantiles <- function(values, stage) {
  probs <- c(0.05, 0.25, 0.5, 0.75, 0.95)
  qs <- apply(values, 2, quantile, probs = probs, na.rm = TRUE)
  tibble(
    stage = stage, sample_id = colnames(values),
    q05 = qs[1, ], q25 = qs[2, ], q50 = qs[3, ], q75 = qs[4, ], q95 = qs[5, ]
  )
}
stopifnot(identical(colnames(before), colnames(dal$data)))
sheets <- list(
  normalized = bind_cols(
    as_tibble(dal$annotation) |> select(uniprot_id, protein, gene, description),
    as_tibble(dal$data)
  ),
  components = components,
  pca_scores = select(pca_scores, sample_id, subject, arm, timepoint, PC1, PC2),
  eta_squared = tibble(uniprot_id = rownames(dal$data), eta2 = unname(eta2)),
  intensity_quantiles = bind_rows(
    intensity_quantiles(before, "before"), intensity_quantiles(as.matrix(dal$data), "after")
  ) |>
    left_join(select(dal$metadata, sample_id, arm, timepoint), by = "sample_id")
)
holds <- c(
  "Cyclic-loess log2 abundance, one row per protein.",
  "Variance explained by the first four principal components.",
  "Each sample's PC1 and PC2 score.",
  "Share of each protein's variance explained by the six group cells.",
  "Each sample's log2 intensity quantiles (5, 25, 50, 75, 95%) before and after cyclic loess."
)
write_workbook(sheets, holds, file.path(out, "normalisation.xlsx"), inputs, paths)

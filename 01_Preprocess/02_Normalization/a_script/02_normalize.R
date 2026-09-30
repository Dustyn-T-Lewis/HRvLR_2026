# Cyclic loess on the filtered DAList, then the blood slope comes off every protein. The matrix
# stays unimputed: limma fits each protein on the samples where it was seen.

suppressPackageStartupMessages({
  library(here)
  library(proteoDA)
  library(dplyr)
  library(tibble)
  library(limma)
  library(ggplot2)
  library(writexl)
})

inputs <- c(filtered = "01_Preprocess/01_Filtering/c_data/DAList_filtered.rds")
paths <- vapply(inputs, here, character(1))
stopifnot(file.exists(paths))
dal <- readRDS(paths[["filtered"]])
out <- here("01_Preprocess", "02_Normalization", "c_data")
reports <- here("01_Preprocess", "02_Normalization", "b_reports")
for (path in c(out, reports)) dir.create(path, recursive = TRUE, showWarnings = FALSE)

# proteoDA draws every method it offers, and QC reports before and after cyclic loess, into
# b_reports/. Cyclic loess is limma's fast method with an adaptive span (chooseLowessSpan), so the
# span moves when filtering changes the protein count.
write_norm_report(dal,
  grouping_column = "group", output_dir = reports,
  filename = "norm_comparison.pdf", overwrite = TRUE
)
write_qc_report(dal,
  color_column = "group", output_dir = reports,
  filename = "qc_pre.pdf", overwrite = TRUE
)
dal <- normalize_data(dal, norm_method = "cycloess")
write_qc_report(dal,
  color_column = "group", output_dir = reports,
  filename = "qc_post.pdf", overwrite = TRUE
)

stopifnot(identical(dal$metadata$sample_id, colnames(dal$data)))
# LR's T3 biopsies carry more blood than any other cell (01_filter.xlsx, blood_by_arm). Each
# protein's slope on the centred blood index is fitted beside the six group means, subject
# blocked, so the slope comes from within-cell spread and no contrast leaks into it. Subtracting
# it leaves the group means as the covariate model estimates them. A protein with no estimable
# slope stays as it was.
dal$metadata$blood <- dal$metadata$blood_rel - mean(dal$metadata$blood_rel)
unadjusted <- dal
blood_design <- model.matrix(~ 0 + group + blood, dal$metadata)
blood_correlation <- duplicateCorrelation(dal$data, blood_design,
  block = dal$metadata$subject
)$consensus.correlation
blood_fit <- lmFit(dal$data, blood_design,
  block = dal$metadata$subject, correlation = blood_correlation
) |>
  eBayes()
slope <- coalesce(blood_fit$coefficients[, "blood"], 0)
dal$data <- dal$data - outer(slope, dal$metadata$blood)
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
saveRDS(unadjusted, file.path(out, "DAList_unadjusted.rds"), compress = "xz")
pdf(file.path(reports, "02_normalize_figures.pdf"), width = 11, height = 8.5)
print(pca_figure)
invisible(dev.off())
sheets <- list(
  normalized = bind_cols(
    as_tibble(dal$annotation) |> select(uniprot_id, protein, gene, description),
    as_tibble(dal$data)
  ),
  components = components,
  pca_scores = select(pca_scores, sample_id, subject, arm, timepoint, PC1, PC2),
  eta_squared = tibble(uniprot_id = rownames(dal$data), eta2 = unname(eta2)),
  blood_slope = tibble(
    uniprot_id = rownames(dal$data), gene = dal$annotation$gene,
    slope = unname(blood_fit$coefficients[, "blood"]), p = unname(blood_fit$p.value[, "blood"])
  ),
  input_manifest = tibble(
    input = names(inputs), path = unname(inputs), md5 = unname(tools::md5sum(paths))
  ),
  package_versions = sessioninfo::package_info("loaded", dependencies = FALSE) |>
    as_tibble() |>
    select(package, version = loadedversion, source)
)
read_me <- tibble(sheet = names(sheets), holds = c(
  "Cyclic-loess log2 abundance with the blood slope removed, one row per protein.",
  "Variance explained by the first four principal components.",
  "Each sample's PC1 and PC2 score.",
  "Share of each protein's variance explained by the six group cells.",
  "log2 change per unit of the centred blood index, with its moderated p, per protein.",
  "Files read, with md5.",
  "Packages loaded at run time."
))
write_xlsx(c(list(read_me = read_me), sheets), file.path(out, "02_normalize.xlsx"))

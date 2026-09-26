# Score every sample on every set with singscore. Scores rest on within-sample ranks, so they do
# not move when the cohort changes. No p-value, no contrast. Ranks need every protein in every
# sample, so this reads the imputed matrix. 05_classify_and_associate_sets consumes the scores.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tibble)
  library(purrr)
  library(ggplot2)
})

out <- here("03_Pathway_Enrichment", "01_Scores", "c_data")
figures <- here("03_Pathway_Enrichment", "01_Scores", "b_reports")
for (path in c(out, figures)) dir.create(path, recursive = TRUE, showWarnings = FALSE)

inputs <- c(
  gene_sets = "03_Pathway_Enrichment/00_Gene_Sets/c_data/gene_sets.rds",
  proteins = "01_Preprocess/03_Imputation/c_data/DAList_imputed.rds"
)
paths <- map_chr(inputs, here)
if (!all(file.exists(paths))) {
  stop(
    "Run 00_Gene_Sets first. Missing: ",
    paste(inputs[!file.exists(paths)], collapse = ", ")
  )
}
gs <- readRDS(paths[["gene_sets"]])
proteins <- readRDS(paths[["proteins"]])

gene_map <- filter(gs$protein_map, selected)
stopifnot(
  identical(gene_map$uniprot_id, intersect(gene_map$uniprot_id, rownames(proteins$data))),
  !anyDuplicated(gene_map$gene)
)

# Sets are keyed on gene symbols.
gene_matrix <- as.matrix(proteins$data)[gene_map$uniprot_id, ]
rownames(gene_matrix) <- gene_map$gene
ranks <- singscore::rankGenes(gene_matrix)
scored <- singscore::multiScore(ranks, upSetColc = gs$sets)
scores <- scored$Scores
stopifnot(
  identical(rownames(scores), names(gs$sets)),
  identical(colnames(scores), colnames(proteins$data))
)
message("scores: ", nrow(scores), " sets x ", ncol(scores), " samples")

# Spread across sets within a sample reflects set composition; spread across samples within a
# set is what the phenotype analysis works with. Both are recorded.
score_summary <- tibble(
  n_sets = nrow(scores), n_samples = ncol(scores),
  min = round(min(scores), 4), max = round(max(scores), 4),
  mean = round(mean(scores), 4),
  mean_sd_across_sets = round(mean(apply(scores, 2, sd)), 4),
  mean_sd_across_samples = round(mean(apply(scores, 1, sd)), 4)
)
print(score_summary)

# Dispersion is the spread of a set's member ranks within one sample: low means the members sit
# together, high that they scatter. multiScore returns it beside the scores at no extra cost.
set_spread <- gs$set_catalog |>
  filter(qualifies) |>
  transmute(set_id, collection) |>
  mutate(
    score = rowMeans(scores[set_id, ]),
    dispersion = rowMeans(scored$Dispersions[set_id, ])
  )
collection_spread <- set_spread |>
  summarise(
    n_sets = n(), median_score = round(median(score), 4),
    median_dispersion = round(median(dispersion), 1), .by = collection
  )
print(as.data.frame(collection_spread))

# How much of the leading components subject identity explains. A large share is why the
# association step reads within-subject change beside the baseline level.
components <- prcomp(t(scores), scale. = FALSE)
variance <- summary(components)$importance[2, 1:4]
targets <- proteins$metadata[match(rownames(components$x), proteins$metadata$sample_id), ]
participant_share <- map_dbl(1:2, function(i) {
  terms <- summary(aov(components$x[, i] ~ targets$subject))[[1]]
  terms[1, 2] / sum(terms[, 2])
})
structure_check <- tibble(
  component = paste0("PC", 1:4),
  variance_explained = round(as.numeric(variance), 3),
  subject_share = c(round(participant_share, 3), NA, NA)
)
print(structure_check)

# singscore's plotDispersion and plotRankDensity draw one signature at a time. These cohort-scale
# views are plain ggplots over the score and dispersion matrices multiScore already returned.
dispersion_figure <- ggplot(set_spread, aes(score, dispersion, colour = collection)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey80") +
  geom_point(alpha = 0.4, size = 0.8) +
  scale_colour_brewer(palette = "Dark2", name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(
    x = "mean score across samples", y = "mean dispersion across samples",
    title = "Set score against dispersion",
    subtitle = sprintf("singscore, %d sets across %d samples", nrow(scores), ncol(scores)),
    caption = paste(
      "One point per set, averaged over samples. Dispersion is the spread of a set's member",
      "ranks within a sample: low means members sit together.",
      "Table: c_data/01_scores.xlsx, set_spread."
    )
  ) +
  theme_minimal(base_size = 9) +
  theme(plot.caption = element_text(size = 7, colour = "grey45", hjust = 0))

group_scores <- tibble(
  group = rep(proteins$metadata$group[match(colnames(scores), proteins$metadata$sample_id)],
    each = nrow(scores)
  ),
  score = as.vector(scores)
)
distribution_figure <- ggplot(group_scores, aes(score, group)) +
  geom_violin(fill = "grey85", colour = NA) +
  geom_boxplot(width = 0.12, outlier.shape = NA, linewidth = 0.3) +
  labs(
    x = "singscore", y = NULL,
    title = "Score distribution by study group",
    subtitle = sprintf("singscore, %d sets pooled", nrow(scores)),
    caption = paste(
      "Every set in every sample, pooled within group. Box is the interquartile range.",
      "Table: c_data/01_scores.xlsx, set_scores."
    )
  ) +
  theme_minimal(base_size = 9) +
  theme(plot.caption = element_text(size = 7, colour = "grey45", hjust = 0))

pdf(file.path(figures, "01_scores_figures.pdf"), width = 11, height = 8.5)
walk(list(dispersion_figure, distribution_figure), print)
invisible(dev.off())

saveRDS(list(scores = scores), file.path(out, "singscore.rds"), compress = "xz")
sheets <- list(
  score_summary = score_summary,
  collection_spread = collection_spread,
  set_spread = set_spread,
  structure_check = structure_check,
  set_scores = rownames_to_column(as.data.frame(scores), "set_id"),
  input_manifest = tibble(
    input = names(inputs), path = unname(inputs), md5 = unname(tools::md5sum(paths))
  ),
  package_versions = sessioninfo::package_info("loaded", dependencies = FALSE) |>
    as_tibble() |>
    select(package, version = loadedversion, source)
)
read_me <- tibble(sheet = names(sheets), holds = c(
  "Range and spread of the score matrix.",
  "Median score and dispersion per collection.",
  "Mean score and dispersion per set across samples.",
  "Variance explained by the first four components, and the subject share of PC1 and PC2.",
  "The set by sample score matrix, as saved in singscore.rds.",
  "Files read, with md5.",
  "Packages loaded at run time."
))
writexl::write_xlsx(c(list(read_me = read_me), sheets), file.path(out, "01_scores.xlsx"))
message("wrote singscore.rds, 01_scores.xlsx and 2 figures")

# Score every sample on every set. singscore ranks each sample's proteins and scores sets
# against those ranks, so a score is rank-based and sample-independent: it does not move when
# the cohort changes, which is what a repeated-measures design needs. No p-value, and it never
# sees the contrast. The matrix is what 05_classify_and_associate_sets reads. Ranks need every
# protein present in every sample, so this reads the imputed matrix.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tibble)
  library(purrr)
  library(ggplot2)
})

out <- here("03_Pathway_Enrichment", "04_run_singscore", "c_data")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

inputs <- c(
  gene_sets = "03_Pathway_Enrichment/00_build_gene_sets/c_data/gene_sets.rds",
  proteins = "01_Preprocess/03_Imputation/c_data/DAList_imputed.rds"
)
paths <- map_chr(inputs, here)
if (!all(file.exists(paths))) {
  stop(
    "Run 00_build_gene_sets first. Missing: ",
    paste(inputs[!file.exists(paths)], collapse = ", ")
  )
}
gs <- readRDS(paths[["gene_sets"]])
proteins <- readRDS(paths[["proteins"]])
manifest <- tibble(
  input = names(inputs), path = unname(inputs), md5 = unname(tools::md5sum(paths))
)

gene_map <- filter(gs$protein_map, selected)
stopifnot(
  identical(gene_map$protein, intersect(gene_map$protein, rownames(proteins$data))),
  !anyDuplicated(gene_map$gene)
)

# Rows become gene symbols because that is the namespace the sets are keyed on.
gene_matrix <- as.matrix(proteins$data)[gene_map$protein, ]
rownames(gene_matrix) <- gene_map$gene
ranks <- singscore::rankGenes(gene_matrix)
scored <- singscore::multiScore(ranks, upSetColc = gs$sets)
scores <- scored$Scores
stopifnot(
  identical(rownames(scores), names(gs$sets)),
  identical(colnames(scores), colnames(proteins$data))
)
message("scores: ", nrow(scores), " sets x ", ncol(scores), " samples")

# Spread across sets within a sample is set composition; spread across samples within a set is
# what the phenotype analysis has to work with, so both are recorded.
score_summary <- tibble(
  sets = nrow(scores), samples = ncol(scores),
  min = round(min(scores), 4), max = round(max(scores), 4),
  mean = round(mean(scores), 4),
  mean_sd_across_sets = round(mean(apply(scores, 2, sd)), 4),
  mean_sd_across_samples = round(mean(apply(scores, 1, sd)), 4)
)
print(score_summary)

# Dispersion is the spread of a set's member ranks within one sample: a low value means the
# members sit together in that sample's ranking, a high one that they are scattered. It comes
# back from multiScore beside the scores, so the cohort view costs no extra call.
set_spread <- gs$set_catalog |>
  filter(qualifies) |>
  transmute(set_id, database) |>
  mutate(
    score = rowMeans(scores[set_id, ]),
    dispersion = rowMeans(scored$Dispersions[set_id, ])
  )
collection_spread <- set_spread |>
  summarise(
    sets = n(), median_score = round(median(score), 4),
    median_dispersion = round(median(dispersion), 1), .by = database
  )
print(as.data.frame(collection_spread))

# How much of the leading components subject identity explains. A large share is why the
# association step reads within-subject change as well as the level.
components <- prcomp(t(scores), scale. = FALSE)
variance <- summary(components)$importance[2, 1:4]
targets <- proteins$metadata[match(rownames(components$x), proteins$metadata$Col_ID), ]
participant_share <- map_dbl(1:2, function(i) {
  terms <- summary(aov(components$x[, i] ~ targets$Subject_ID))[[1]]
  terms[1, 2] / sum(terms[, 2])
})
structure_check <- tibble(
  component = paste0("PC", 1:4),
  variance_explained = round(as.numeric(variance), 3),
  subject_share = c(round(participant_share, 3), NA, NA)
)
print(structure_check)

# singscore ships plotDispersion and plotRankDensity, which draw one signature at a time.
# Reporting here is at cohort scale, so these are plain ggplots over the two matrices multiScore
# already returned.
figures <- here("03_Pathway_Enrichment", "04_run_singscore", "b_reports")
dir.create(figures, recursive = TRUE, showWarnings = FALSE)
save_figure <- function(figure, name, height) {
  walk(c("png", "pdf"), \(extension) {
    ggsave(file.path(figures, paste0(name, ".", extension)), figure,
      width = 6.5, height = height, dpi = 200, bg = "white"
    )
  })
}
save_figure(
  ggplot(set_spread, aes(score, dispersion, colour = database)) +
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
        "ranks within a sample: low means members sit together. Table: c_data/set_scores.csv."
      )
    ) +
    theme_minimal(base_size = 9) +
    theme(plot.caption = element_text(size = 7, colour = "grey45", hjust = 0)),
  "01_score_dispersion",
  height = 4.5
)

group_scores <- tibble(
  group = rep(proteins$metadata$Group_Time[match(colnames(scores), proteins$metadata$Col_ID)],
    each = nrow(scores)
  ),
  score = as.vector(scores)
)
save_figure(
  ggplot(group_scores, aes(score, group)) +
    geom_violin(fill = "grey85", colour = NA) +
    geom_boxplot(width = 0.12, outlier.shape = NA, linewidth = 0.3) +
    labs(
      x = "singscore", y = NULL,
      title = "Score distribution by study group",
      subtitle = sprintf("singscore, %d sets pooled", nrow(scores)),
      caption = paste(
        "Every set in every sample, pooled within group. Box is the interquartile range.",
        "Scores are rank-based and sample-independent. Table: c_data/set_scores.csv."
      )
    ) +
    theme_minimal(base_size = 9) +
    theme(plot.caption = element_text(size = 7, colour = "grey45", hjust = 0)),
  "02_score_distribution",
  height = 3.5
)
message("wrote 2 cohort figures")

packages <- c("here", "proteoDA", "singscore", "dplyr", "purrr", "ggplot2")
versions <- tibble(
  package = packages, version = map_chr(packages, \(p) as.character(packageVersion(p)))
)
score_table <- rownames_to_column(as.data.frame(scores), "set_id")

saveRDS(
  list(
    scores = scores, score_summary = score_summary, structure_check = structure_check,
    collection_spread = collection_spread,
    provenance = list(
      created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
      inputs = manifest, packages = versions
    )
  ),
  file.path(out, "singscore.rds"),
  compress = "xz"
)
writexl::write_xlsx(
  list(
    score_summary = score_summary,
    collection_spread = collection_spread,
    structure_check = structure_check,
    set_scores = score_table,
    input_manifest = manifest,
    package_versions = versions
  ),
  file.path(out, "04_run_singscore.xlsx")
)
readr::write_csv(score_table, file.path(out, "set_scores.csv"))
combined <- file.path(figures, "04_run_singscore_figures.pdf")
pages <- setdiff(list.files(figures, "[.]pdf$", full.names = TRUE), combined)
invisible(qpdf::pdf_combine(sort(pages), combined))
message(
  "wrote singscore.rds, 04_run_singscore.xlsx, set_scores.csv and a ",
  length(pages), "-page figure PDF"
)

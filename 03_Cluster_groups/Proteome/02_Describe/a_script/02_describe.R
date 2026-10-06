# What separates each proteome split, at the input that made it. Per input: each module
# eigengene's mean (z-scored) per cluster, limma on every protein at that input (cluster 2 minus
# cluster 1), and over-representation (clusterProfiler::enricher) among the proteins higher and
# lower in cluster 2 at p < 0.05. The clusters were built from these data, so every p here is
# descriptive: it says what the split looks like, not that the split is real.

suppressPackageStartupMessages({
  library(here)
  library(limma)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
})
source(here("R", "helpers.R"))

stage <- here("03_Cluster_groups", "Proteome", "02_Describe")
out <- file.path(stage, "c_data")
figures <- file.path(stage, "b_reports")
walk(c(out, figures), dir.create, recursive = TRUE, showWarnings = FALSE)
inputs <- c(
  clusters = "03_Cluster_groups/Proteome/01_Cluster/c_data/01_cluster.xlsx",
  corrected = "01_Trunk/03_Quantify/c_data/DAList_corrected.rds",
  modules = "01_Trunk/06_Modules/c_data/modules.rds",
  gene_sets = "01_Trunk/04_Gene_Sets/c_data/gene_sets.rds"
)
paths <- map_chr(inputs, here)
stopifnot(file.exists(paths))

assignments <- readxl::read_excel(paths[["clusters"]], "assignments")
dal <- readRDS(paths[["corrected"]])
meta <- as_tibble(dal$metadata)
modules <- readRDS(paths[["modules"]])
gs <- readRDS(paths[["gene_sets"]])
gene_of <- with(filter(gs$protein_map, selected), set_names(gene, uniprot_id))
protein_label <- with(dal$annotation, set_names(gene, uniprot_id))
set_label <- with(gs$set_catalog, set_names(pathway, set_id))
term2gene <- tibble(set_id = rep(names(gs$sets), lengths(gs$sets)), gene = unlist(gs$sets))
input_matrix <- function(values, input) {
  switch(input,
    T1 = subject_values(values, meta, "T1"),
    T2 = subject_values(values, meta, "T2"),
    T3 = subject_values(values, meta, "T3"),
    training = subject_change(values, meta, "T1", "T2"),
    acute = subject_change(values, meta, "T2", "T3")
  )
}
cluster_of <- \(input) with(filter(assignments, input == .env$input), set_names(cluster, subject))

module_means <- map(unique(assignments$input), \(input) {
  x <- input_matrix(modules$eigengenes, input)
  z <- t(scale(t(x)))
  clusters <- cluster_of(input)[colnames(z)]
  as_tibble(z, rownames = "module") |>
    pivot_longer(-module, names_to = "subject", values_to = "z") |>
    mutate(cluster = clusters[subject]) |>
    summarise(mean_z = mean(z), .by = c(module, cluster)) |>
    mutate(input = input, .before = 1)
}) |>
  list_rbind()

protein_difference <- map(unique(assignments$input), \(input) {
  x <- input_matrix(dal$data, input)
  predictors <- data.frame(cluster = factor(cluster_of(input)[colnames(x)]))
  subject_fit(x, predictors, ~cluster, "cluster2") |>
    mutate(input = input, label = unname(protein_label[feature]), .before = 1) |>
    mutate(fdr = p.adjust(p, "BH"))
}) |>
  list_rbind()

ora <- function(members, universe) {
  res <- clusterProfiler::enricher(na.omit(unname(gene_of[members])),
    universe = na.omit(unname(gene_of[universe])), TERM2GENE = term2gene, minGSSize = 15,
    maxGSSize = 500, pvalueCutoff = 1, qvalueCutoff = 1
  )
  if (is.null(res)) {
    return(tibble())
  }
  as_tibble(res@result) |>
    transmute(
      set_id = ID, pathway = unname(set_label[ID]), overlap = Count, p = pvalue,
      fdr = p.adjust
    )
}
enrichment <- protein_difference |>
  filter(p < 0.05) |>
  mutate(direction = if_else(effect > 0, "higher in cluster 2", "lower in cluster 2")) |>
  nest(.by = c(input, direction)) |>
  mutate(result = map(data, \(d) ora(d$feature, rownames(dal$data)))) |>
  select(-data) |>
  unnest(result) |>
  arrange(input, direction, p)

pdf(file.path(figures, "02_describe_figures.pdf"), width = 11, height = 8.5)
print(
  ggplot(module_means, aes(factor(cluster), module, fill = mean_z)) +
    geom_tile() +
    facet_wrap(~input, nrow = 1) +
    scale_fill_gradient2(low = "#2166AC", high = "#B2182B", name = "mean z") +
    labs(
      title = "Module eigengenes per proteome cluster, at the input that made the split",
      subtitle = "Descriptive: the clusters were built from these values", x = "Cluster", y = NULL
    ) +
    theme_bw(base_size = 9)
)
print(
  enrichment |>
    slice_min(p, n = 8, by = c(input, direction), with_ties = FALSE) |>
    ggplot(aes(-log10(p), reorder(stringr::str_trunc(pathway, 50), -log10(p)), fill = direction)) +
    geom_col() +
    facet_wrap(~input, scales = "free_y", ncol = 3) +
    labs(
      title = "Top over-represented sets among proteins differing between clusters (p < 0.05)",
      subtitle = "Descriptive; enricher against all measured proteins", x = "-log10 p", y = NULL
    ) +
    theme_bw(base_size = 7) +
    theme(legend.position = "bottom")
)
invisible(dev.off())

write_workbook(
  list(
    module_means = module_means, protein_difference = protein_difference, enrichment = enrichment
  ),
  c(
    "Per input, module and cluster: mean z-scored eigengene.",
    "Per input and protein: limma cluster 2 minus cluster 1, t, p, BH (descriptive).",
    "Per input and direction: enricher sets among p < 0.05 proteins, overlap, p, BH (descriptive)."
  ),
  file.path(out, "02_describe.xlsx"), inputs, paths
)

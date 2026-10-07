# Proteins grouped by how they relate to phenotype, every subject pooled, no group label.
#   profiles   each protein described by its limma t against every trait in every pairing
#              (01_Trait_Associations/01_Proteins); average linkage on 1 - Pearson r, cut by
# dynamicTreeCut
#   per trait  for each trait's training change and T1 level, the proteins at p < 0.05 grouped by
#              how they co-vary across people, with each group's direction and main module
# Every group is described by over-representation (clusterProfiler::enricher) against the proteins
# it was drawn from. The groups come from the same t values they summarise, so the descriptions are
# not tests of the grouping; the mean silhouette reports how separate the groups are.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(cluster)
})
source(here("R", "helpers.R"))

stage <- here("04_Continuous_no_group", "03_Profile_Clusters")
out <- file.path(stage, "c_data")
figures <- file.path(stage, "b_reports")
walk(c(out, figures), dir.create, recursive = TRUE, showWarnings = FALSE)
inputs <- c(
  corrected = "01_Trunk/03_Quantification/c_data/DAList_corrected.rds",
  associate = "04_Continuous_no_group/01_Trait_Associations/01_Proteins/c_data/associate.rds",
  gene_sets = "01_Trunk/04_Gene_Sets/c_data/gene_sets.rds",
  modules = "01_Trunk/06_Module_Building/c_data/modules.rds"
)
paths <- map_chr(inputs, here)
stopifnot(file.exists(paths))

dal <- readRDS(paths[["corrected"]])
meta <- as_tibble(dal$metadata)
gs <- readRDS(paths[["gene_sets"]])
modules <- readRDS(paths[["modules"]])
gene_of <- with(filter(gs$protein_map, selected), set_names(gene, uniprot_id))
protein_label <- with(dal$annotation, set_names(gene, uniprot_id))
set_label <- with(gs$set_catalog, set_names(pathway, set_id))
module_of <- with(modules$membership, set_names(module, uniprot_id))
observed_t <- readRDS(paths[["associate"]])$change_score |> select(window, trait, feature, t, p)
term2gene <- tibble(set_id = rep(names(gs$sets), lengths(gs$sets)), gene = unlist(gs$sets))

tree_groups <- function(m, min_size) {
  d <- as.dist(1 - cor(t(m)))
  groups <- dynamicTreeCut::cutreeDynamic(hclust(d, "average"),
    distM = as.matrix(d), minClusterSize = min_size, deepSplit = 2, verbose = 0
  )
  kept <- groups != 0
  silhouette_mean <- if (n_distinct(groups[kept]) < 2) {
    NA_real_
  } else {
    mean(silhouette(groups[kept], as.dist(as.matrix(d)[kept, kept]))[, "sil_width"])
  }
  list(groups = set_names(groups, rownames(m))[kept], silhouette = silhouette_mean)
}
describe <- function(members, universe) {
  res <- clusterProfiler::enricher(na.omit(unname(gene_of[members])),
    universe = na.omit(unname(gene_of[universe])), TERM2GENE = term2gene, minGSSize = 15,
    maxGSSize = 500, pvalueCutoff = 1, qvalueCutoff = 1
  )
  if (is.null(res) || nrow(res@result) == 0) {
    return(tibble(top_set = NA_character_, top_set_p = NA_real_, top_set_fdr = NA_real_))
  }
  top <- slice_min(as_tibble(res@result), pvalue, n = 1, with_ties = FALSE)
  tibble(top_set = unname(set_label[top$ID]), top_set_p = top$pvalue, top_set_fdr = top$p.adjust)
}
summarise_groups <- function(groups, t, universe) {
  tibble(feature = names(groups), group = groups, t = t[names(groups)]) |>
    summarise(
      n_proteins = n(), mean_t = mean(t), share_positive = mean(t > 0), members = list(feature),
      main_module = names(which.max(table(module_of[feature]))),
      proteins = paste(head(protein_label[feature][order(-abs(t))], 15), collapse = ", "),
      .by = group
    ) |>
    mutate(enrichment = map(members, \(m) describe(m, universe))) |>
    unnest(enrichment) |>
    select(-members) |>
    arrange(group)
}

profiles <- observed_t |>
  transmute(feature, screen = paste(window, trait), t) |>
  pivot_wider(names_from = screen, values_from = t) |>
  column_to_rownames("feature") |>
  as.matrix() |>
  na.omit()
profile_cut <- tree_groups(profiles, min_size = 20)
profile_summary <- summarise_groups(
  profile_cut$groups, rowMeans(profiles), rownames(profiles)
)
profile_means <- as_tibble(profiles[names(profile_cut$groups), ], rownames = "feature") |>
  mutate(group = profile_cut$groups) |>
  pivot_longer(-c(feature, group), names_to = "screen", values_to = "t") |>
  summarise(mean_t = mean(t), .by = c(group, screen))

matrix_of <- list(
  training = subject_change(dal$data, meta, "T1", "T2"),
  baseline_level = subject_values(dal$data, meta, "T1")
)
trait_groups <- observed_t |>
  filter(window %in% names(matrix_of), p < 0.05) |>
  nest(.by = c(window, trait)) |>
  mutate(result = map2(window, data, \(w, hits) {
    cut <- tree_groups(matrix_of[[w]][hits$feature, , drop = FALSE], min_size = 8)
    universe <- filter(observed_t, window == w)$feature
    summarise_groups(cut$groups, set_names(hits$t, hits$feature), universe) |>
      mutate(silhouette = cut$silhouette, n_hits = nrow(hits))
  })) |>
  select(-data) |>
  unnest(result)

pdf(file.path(figures, "profile_clusters_report.pdf"), width = 11, height = 8.5)
print(
  ggplot(profile_means, aes(screen, factor(group), fill = mean_t)) +
    geom_tile() +
    scale_fill_gradient2(low = "#2166AC", high = "#B2182B", name = "mean t") +
    labs(
      title = "Protein association profiles: mean limma t per group and screen",
      subtitle = sprintf(
        "%d proteins in %d groups; mean silhouette %.2f",
        length(profile_cut$groups), n_distinct(profile_cut$groups), profile_cut$silhouette
      ),
      x = NULL, y = "Profile group"
    ) +
    theme_bw(base_size = 7) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
)
print(
  ggplot(trait_groups, aes(n_proteins, paste(window, trait), fill = mean_t)) +
    geom_col(position = position_stack()) +
    scale_fill_gradient2(low = "#2166AC", high = "#B2182B", name = "mean t") +
    labs(
      title = "Co-variation groups among each trait's p < 0.05 proteins",
      subtitle = "One bar segment per group; fill: the group's mean t against the trait",
      x = "Proteins", y = NULL
    ) +
    theme_bw(base_size = 8)
)
invisible(dev.off())

write_workbook(
  list(
    profile_groups = profile_summary,
    profile_members = tibble(
      feature = names(profile_cut$groups), label = unname(protein_label[names(profile_cut$groups)]),
      group = profile_cut$groups
    ),
    profile_means = profile_means,
    trait_groups = trait_groups
  ),
  c(
    sprintf(
      "Profile groups: size, mean t, main module, top over-represented set; mean silhouette %.3f.",
      profile_cut$silhouette
    ),
    "Each protein's profile group.",
    "Mean limma t per profile group and screen.",
    "Per trait and window: co-variation groups among p < 0.05 proteins, with silhouette and ORA."
  ),
  file.path(out, "profile_clusters.xlsx"), inputs, paths
)

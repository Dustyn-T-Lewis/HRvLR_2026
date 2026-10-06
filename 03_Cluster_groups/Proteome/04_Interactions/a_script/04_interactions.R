# Two interactions with the proteome clusters, at protein and module level.
#   cluster x time       clusters made at one input, tested on a change that comes later and was
#                        not used to make them: T1 clusters on training and acute change, T2 and
#                        training clusters on acute change. limma, change ~ cluster.
#   cluster x phenotype  does a feature's training change track a trait differently in the two
#                        clusters? limma, change ~ cluster * trait, the interaction term, per trait.
# The three-way cluster x time x phenotype model is not fitted: with 14 or 15 people it leaves too
# few residual degrees of freedom to estimate. BH within each screen.

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

stage <- here("03_Cluster_groups", "Proteome", "04_Interactions")
out <- file.path(stage, "c_data")
figures <- file.path(stage, "b_reports")
walk(c(out, figures), dir.create, recursive = TRUE, showWarnings = FALSE)
inputs <- c(
  clusters = "03_Cluster_groups/Proteome/01_Cluster/c_data/01_cluster.xlsx",
  corrected = "01_Trunk/03_Quantify/c_data/DAList_corrected.rds",
  modules = "01_Trunk/06_Modules/c_data/modules.rds",
  phenotype = "00_Input/phenotype.csv"
)
paths <- map_chr(inputs, here)
stopifnot(file.exists(paths))

assignments <- readxl::read_excel(paths[["clusters"]], "assignments")
dal <- readRDS(paths[["corrected"]])
meta <- as_tibble(dal$metadata)
modules <- readRDS(paths[["modules"]])
phenotype <- readr::read_csv(paths[["phenotype"]], show_col_types = FALSE)
training_traits <- phenotype_tables(phenotype)$pairings$training
cluster_of <- \(input) with(filter(assignments, input == .env$input), set_names(cluster, subject))
label_of <- c(
  with(dal$annotation, set_names(gene, uniprot_id)), set_names(rownames(modules$eigengenes))
)
changes <- list(
  proteins = list(
    training = subject_change(dal$data, meta, "T1", "T2"),
    acute = subject_change(dal$data, meta, "T2", "T3")
  ),
  modules = list(
    training = subject_change(modules$eigengenes, meta, "T1", "T2"),
    acute = subject_change(modules$eigengenes, meta, "T2", "T3")
  )
)
later <- tribble(
  ~input, ~window,
  "T1", "training", "T1", "acute", "T2", "acute", "training", "acute"
)

cluster_time <- expand_grid(later, level = names(changes)) |>
  mutate(result = pmap(list(input, window, level), \(input, window, level) {
    x <- changes[[level]][[window]]
    clusters <- cluster_of(input)[colnames(x)]
    keep <- !is.na(clusters)
    subject_fit(
      x[, keep, drop = FALSE], data.frame(cluster = factor(clusters[keep])), ~cluster, "cluster2"
    ) |>
      mutate(fdr = p.adjust(p, "BH"))
  })) |>
  unnest(result) |>
  mutate(label = unname(label_of[feature]), .after = feature)

cluster_phenotype <- expand_grid(
  input = unique(assignments$input), level = names(changes),
  trait = setdiff(names(training_traits), "subject")
) |>
  mutate(result = pmap(list(input, level, trait), \(input, level, trait) {
    x <- changes[[level]][["training"]]
    predictors <- data.frame(
      cluster = factor(cluster_of(input)[colnames(x)]),
      trait = pairing_value(list(training = training_traits), "training", trait, colnames(x))
    )
    subject_fit(x, predictors, ~ cluster * trait, "cluster2:trait") |>
      mutate(fdr = p.adjust(p, "BH"))
  })) |>
  unnest(result) |>
  mutate(label = unname(label_of[feature]), .after = feature)

counts <- bind_rows(
  summarise(cluster_time,
    n_tested = n(), n_nominal = sum(p < 0.05), n_fdr_05 = sum(fdr < 0.05),
    .by = c(input, window, level)
  ) |>
    mutate(analysis = "cluster x time", trait = NA_character_),
  summarise(cluster_phenotype,
    n_tested = n(), n_nominal = sum(p < 0.05), n_fdr_05 = sum(fdr < 0.05),
    .by = c(input, level, trait)
  ) |>
    mutate(analysis = "cluster x phenotype", window = "training")
) |>
  relocate(analysis, input, window, trait, level)

pdf(file.path(figures, "04_interactions_figures.pdf"), width = 11, height = 8.5)
print(
  ggplot(counts, aes(n_nominal / n_tested, paste(input, "->", coalesce(trait, window)))) +
    geom_col(aes(fill = level), position = position_dodge(width = 0.8), width = 0.75) +
    geom_vline(xintercept = 0.05, linetype = "dashed", linewidth = 0.3) +
    facet_wrap(~analysis, scales = "free_y") +
    scale_x_continuous(labels = scales::label_percent()) +
    labs(
      title = "Share of features at p < 0.05 per interaction screen",
      subtitle = "Dashed: 5%. Cluster input -> later window, or -> trait for cluster x phenotype",
      x = "Features at p < 0.05", y = NULL, fill = NULL
    ) +
    theme_bw(base_size = 8)
)
invisible(dev.off())

write_workbook(
  list(
    counts = counts,
    cluster_time = arrange(cluster_time, input, window, level, p),
    cluster_phenotype = arrange(cluster_phenotype, input, level, trait, p)
  ),
  c(
    "Per screen: features tested, at p < 0.05 and at BH < 0.05.",
    "Clusters from an earlier input on a later change: cluster 2 minus 1, t, p, BH.",
    "Training change ~ cluster * trait: the interaction (cluster 2 slope minus 1), t, p, BH."
  ),
  file.path(out, "04_interactions.xlsx"), inputs, paths
)

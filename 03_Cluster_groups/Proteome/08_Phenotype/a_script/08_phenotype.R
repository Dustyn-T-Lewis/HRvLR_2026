# Do the proteome clusters differ in size or strength? Phenotype played no part in the clustering,
# so these tests are not circular. Per input and trait: each cluster's median, the Wilcoxon
# rank-sum test (wilcox.test, exact where it can be) and BH within input. Traits: % change over
# training, the composite and volume load, and the T1 and T2 levels. The cross-table with HR/LR is
# a description, with Fisher's exact p.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
})
source(here("R", "helpers.R"))

stage <- here("03_Cluster_groups", "Proteome", "03_Phenotype")
out <- file.path(stage, "c_data")
figures <- file.path(stage, "b_reports")
walk(c(out, figures), dir.create, recursive = TRUE, showWarnings = FALSE)
inputs <- c(
  clusters = "03_Cluster_groups/Proteome/01_Cluster/c_data/01_cluster.xlsx",
  phenotype = "00_Input/phenotype.csv"
)
paths <- map_chr(inputs, here)
stopifnot(file.exists(paths))

assignments <- readxl::read_excel(paths[["clusters"]], "assignments")
phenotype <- readr::read_csv(paths[["phenotype"]], show_col_types = FALSE)
pairings <- phenotype_tables(phenotype)$pairings
traits <- bind_rows(
  pivot_longer(pairings$training, -subject, names_to = "trait") |> mutate(kind = "change"),
  pivot_longer(pairings$baseline_level, -subject, names_to = "trait") |> mutate(kind = "T1 level"),
  pivot_longer(pairings$post_level, -subject, names_to = "trait") |> mutate(kind = "T2 level")
)

tests <- assignments |>
  select(input, subject, cluster) |>
  inner_join(traits, by = "subject", relationship = "many-to-many") |>
  filter(!is.na(value)) |>
  nest(.by = c(input, kind, trait)) |>
  mutate(test = map(data, \(d) {
    rank_sum <- suppressWarnings(wilcox.test(value ~ cluster, data = d))
    tibble(
      n_1 = sum(d$cluster == 1), n_2 = sum(d$cluster == 2),
      median_1 = median(d$value[d$cluster == 1]), median_2 = median(d$value[d$cluster == 2]),
      w = unname(rank_sum$statistic), p = rank_sum$p.value
    )
  })) |>
  select(-data) |>
  unnest(test) |>
  mutate(fdr = p.adjust(p, "BH"), .by = input) |>
  arrange(input, p)

hr_lr <- assignments |>
  summarise(
    hr_in_1 = sum(arm == "HR" & cluster == 1), lr_in_1 = sum(arm == "LR" & cluster == 1),
    hr_in_2 = sum(arm == "HR" & cluster == 2), lr_in_2 = sum(arm == "LR" & cluster == 2),
    fisher_p = fisher.test(table(arm, cluster))$p.value,
    .by = input
  )

pdf(file.path(figures, "03_phenotype_figures.pdf"), width = 11, height = 8.5)
print(
  assignments |>
    inner_join(filter(traits, kind == "change"), by = "subject", relationship = "many-to-many") |>
    ggplot(aes(factor(cluster), value)) +
    geom_boxplot(outlier.shape = NA, width = 0.5) +
    geom_jitter(width = 0.1, size = 1) +
    facet_grid(trait ~ input, scales = "free_y") +
    labs(
      title = "Training response by proteome cluster",
      subtitle = "Wilcoxon p per panel in c_data/03_phenotype.xlsx, tests",
      x = "Cluster", y = "Value"
    ) +
    theme_bw(base_size = 7)
)
invisible(dev.off())

write_workbook(
  list(tests = tests, hr_lr = hr_lr),
  c(
    "Per input, trait kind and trait: n and median per cluster, Wilcoxon W and p, BH within input.",
    "Per input: HR and LR counts in each cluster, Fisher's exact p (description only)."
  ),
  file.path(out, "03_phenotype.xlsx"), inputs, paths
)

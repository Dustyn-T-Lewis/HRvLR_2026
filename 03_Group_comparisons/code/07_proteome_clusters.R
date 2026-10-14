#' ---
#' title: "Proteome clusters: phenotype, and every result stacked"
#' params:
#'   grouping: Proteome_clusters
#' format:
#'   typst:
#'     toc: true
#'     margin:
#'       x: 1.6cm
#' execute:
#'   echo: true
#'   warning: false
#' ---
#'
#' Do the people the proteome clustered together differ in muscle size or strength? The phenotype
#' played no part in the clustering, so these tests are not circular. This step reads all five
#' proteome groupings, so it runs after their comparisons, and ends by stacking their results.
#'
#' For each clustering input and trait: each cluster's median, the Wilcoxon rank-sum test (exact
#' where `wilcox.test()` can compute it) and BH within input. Traits are the % change over training,
#' the composite and volume load, and the T1 and T2 levels. The cross-table with HR/LR describes,
#' with Fisher's exact p; it decides nothing.
#'
#' ## Packages

#| label: packages
#| message: false
library(here)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(ggplot2)
source(here("00_Setup", "Shared_functions", "analysis_reporting_helpers.R"))
source(here("00_Setup", "Shared_functions", "proteome_inputs.R"))

#' ## Inputs

#| label: inputs
inputs <- c(
  set_names(
    file.path("02_Groups", paste0("Proteome_", inputs_used), "c_data/clustering.xlsx"),
    paste0("clusters_", inputs_used)
  ),
  phenotype = "00_Setup/Input/phenotype.csv"
)
paths <- map_chr(inputs, here)
stopifnot(file.exists(paths))
assignments <- map(paste0("clusters_", inputs_used), \(name) {
  readxl::read_excel(paths[[name]], "assignments")
}) |>
  list_rbind()
phenotype <- readr::read_csv(paths[["phenotype"]], show_col_types = FALSE)
pairings <- phenotype_tables(phenotype)$pairings
traits <- bind_rows(
  pivot_longer(pairings$training, -subject, names_to = "trait") |> mutate(kind = "change"),
  pivot_longer(pairings$baseline_level, -subject, names_to = "trait") |> mutate(kind = "T1 level"),
  pivot_longer(pairings$post_level, -subject, names_to = "trait") |> mutate(kind = "T2 level")
)
count(assignments, input, cluster)

#' ## Tests

#| label: tests
tests <- assignments |>
  select(input, subject, cluster) |>
  inner_join(traits, by = "subject", relationship = "many-to-many") |>
  filter(!is.na(value)) |>
  nest(.by = c(input, kind, trait)) |>
  mutate(test = map(data, \(x) {
    rank_sum <- suppressWarnings(wilcox.test(value ~ cluster, data = x))
    tibble(
      n_1 = sum(x$cluster == 1), n_2 = sum(x$cluster == 2),
      median_1 = median(x$value[x$cluster == 1]), median_2 = median(x$value[x$cluster == 2]),
      w = unname(rank_sum$statistic), p = rank_sum$p.value
    )
  })) |>
  select(-data) |>
  unnest(test) |>
  mutate(fdr = p.adjust(p, "BH"), .by = input) |>
  arrange(input, p)
head(tests)


#| label: hr-lr
hr_lr <- assignments |>
  summarise(
    hr_in_1 = sum(arm == "HR" & cluster == 1), lr_in_1 = sum(arm == "LR" & cluster == 1),
    hr_in_2 = sum(arm == "HR" & cluster == 2), lr_in_2 = sum(arm == "LR" & cluster == 2),
    fisher_p = fisher.test(table(arm, cluster))$p.value,
    .by = input
  )
hr_lr


#| label: fig-response
#| fig-cap: "Training response by proteome cluster, per input."
#| fig-width: 11
#| fig-height: 8.5
assignments |>
  inner_join(filter(traits, kind == "change"), by = "subject", relationship = "many-to-many") |>
  ggplot(aes(factor(cluster), value)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_point(position = position_jitter(width = 0.1, seed = 1), size = 1) +
  facet_grid(trait ~ input, scales = "free_y") +
  labs(x = "cluster", y = "value") +
  theme_bw(base_size = 7)

#' ## Validation

#| label: validate
stopifnot(all(tests$fdr >= tests$p), all(tests$n_1 >= 1), all(tests$n_2 >= 1))

#' ## Summary
#'
#' `r nrow(tests)` cluster-against-trait tests over `r n_distinct(tests$input)` inputs:
#' `r sum(tests$p < 0.05)` at p < 0.05, `r sum(tests$fdr < 0.05)` at BH < 0.05.
#'
#' ## Save

#| label: save
out <- here("03_Group_comparisons", params$grouping, "c_data")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write_workbook(
  list(tests = tests, hr_lr = hr_lr),
  c(
    "Per input, trait kind and trait: n and median per cluster, Wilcoxon W and p, BH within input.",
    "Per input: HR and LR counts in each cluster, Fisher's exact p (description only)."
  ),
  file.path(out, "07_cluster_phenotype.xlsx"), inputs, paths
)

#' ## Every result, stacked
#'
#' Each proteome grouping writes its own workbooks. `stack_inputs()` (`proteome_inputs.R`) stacks
#' one workbook's result sheets from all five into `<file>_all_inputs.xlsx`, with an `input` column
#' and a `validity` column: a contrast that uses exactly the time points the clustering was built on
#' is `circular`, one that shares a time point is `caution`, the rest `clean`. The trait
#' associations workbook leaves out its within-group sheet, which would pass the xlsx row limit.

#| label: stack
walk(
  c(
    "01_design_contrasts.xlsx", "02_protein_contrasts.xlsx", "03_set_contrasts.xlsx",
    "04_module_contrasts.xlsx", "05_protein_classification.xlsx", "05_set_classification.xlsx",
    "05_module_classification.xlsx"
  ),
  stack_inputs
)
stack_inputs("06_trait_associations.xlsx", leave_out = "within_group")

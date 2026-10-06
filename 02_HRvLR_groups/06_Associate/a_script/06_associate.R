# Phenotype links inside each group and between groups, for proteins, gene-set scores and module
# eigengenes. The pooled, label-free associations live in 04_Continuous_no_group; this stage asks
# only what the grouping adds.
#
# Within group: Spearman between each feature and each outcome, exact p (cor.test), in each group.
# Group difference: limma on feature ~ group * outcome, the second group's slope minus the first's.
# Flips: DGCA's Fisher z difference between the groups' Spearman correlations, with DGCA's own p
# for that difference (pValDiff).
# The outcomes that define the grouping stay out of the group difference and flips: they re-measure
# the split. Each group needs five subjects with the outcome measured.

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

stage <- here("02_HRvLR_groups", "06_Associate")
out <- file.path(stage, "c_data")
figures <- file.path(stage, "b_reports")
walk(c(out, figures), dir.create, recursive = TRUE, showWarnings = FALSE)

inputs <- c(
  corrected = "01_Trunk/03_Quantify/c_data/DAList_corrected.rds",
  singscore = "01_Trunk/05_Scores/c_data/singscore.rds",
  gene_sets = "01_Trunk/04_Gene_Sets/c_data/gene_sets.rds",
  modules = "01_Trunk/06_Modules/c_data/modules.rds",
  phenotype = "00_Input/phenotype.csv"
)
paths <- map_chr(inputs, here)
stopifnot(file.exists(paths))

phenotype <- readr::read_csv(paths[["phenotype"]], show_col_types = FALSE)
group_of <- set_names(phenotype$arm, phenotype$subject)
group_levels <- c("LR", "HR")
defining_outcomes <- "comp_hypertrophy"

dal <- readRDS(paths[["corrected"]])
catalog <- readRDS(paths[["gene_sets"]])$set_catalog |> filter(qualifies)
modules <- readRDS(paths[["modules"]])
meta <- as_tibble(dal$metadata)
feature_levels <- list(
  proteins = list(values = dal$data, label = with(dal$annotation, set_names(gene, uniprot_id))),
  sets = list(
    values = readRDS(paths[["singscore"]])$scores[catalog$set_id, ],
    label = with(catalog, set_names(pathway, set_id)),
    collection = with(catalog, set_names(collection, set_id))
  ),
  modules = list(
    values = modules$eigengenes, label = set_names(rownames(modules$eigengenes))
  )
)
walk(feature_levels, \(level) stopifnot(identical(colnames(level$values), meta$sample_id)))

pairings <- phenotype_tables(phenotype)$pairings
outcomes_of <- map(pairings, \(tt) setdiff(names(tt), "subject"))
window_of <- \(values) {
  list(
    training = subject_change(values, meta, "T1", "T2"),
    baseline = subject_values(values, meta, "T1"),
    acute = subject_change(values, meta, "T2", "T3"),
    baseline_level = subject_values(values, meta, "T1"),
    post_level = subject_values(values, meta, "T2"),
    post_bout_level = subject_values(values, meta, "T3")
  )
}
enough <- \(subjects, outcome) min(table(group_of[subjects][!is.na(outcome)])) >= 5

per_level <- imap(feature_levels, \(feature_level, level_name) {
  windows <- window_of(feature_level$values)
  screens <- imap(windows, \(x, window) {
    map(set_names(outcomes_of[[window]]), \(name) {
      list(x = x, outcome = pairing_value(pairings, window, name, colnames(x)))
    })
  })

  within_group <- imap(screens, \(by_outcome, window) {
    imap(by_outcome, \(s, name) {
      map(set_names(group_levels), \(g) {
        kept <- group_of[colnames(s$x)] == g
        spearman_by_row(s$x[, kept, drop = FALSE], s$outcome[kept], exact = TRUE)
      }) |>
        list_rbind(names_to = "group") |>
        mutate(outcome = name)
    }) |>
      list_rbind() |>
      mutate(window = window)
  }) |>
    list_rbind()

  group_difference <- imap(screens, \(by_outcome, window) {
    by_outcome[setdiff(names(by_outcome), defining_outcomes)] |>
      keep(\(s) enough(colnames(s$x), s$outcome)) |>
      imap(\(s, name) {
        predictors <- data.frame(
          group = factor(group_of[colnames(s$x)], group_levels), outcome = s$outcome
        )
        term <- sprintf("group%s:outcome", group_levels[2])
        subject_fit(s$x, predictors, ~ group * outcome, term) |>
          mutate(outcome = name)
      }) |>
      list_rbind() |>
      mutate(window = window)
  }) |>
    list_rbind() |>
    select(-coefficient)

  flips <- imap(screens, \(by_outcome, window) {
    by_outcome[setdiff(names(by_outcome), defining_outcomes)] |>
      keep(\(s) enough(colnames(s$x), s$outcome)) |>
      imap(\(s, name) {
        ok <- !is.na(s$outcome)
        group <- factor(group_of[colnames(s$x)][ok], rev(group_levels))
        trait_row <- matrix(s$outcome[ok], nrow = 1, dimnames = list(name, NULL))
        design <- DGCA::makeDesign(group)
        observed <- DGCA::ddcorAll(s$x[, ok, drop = FALSE], design, levels(group),
          inputMatB = trait_row, corrType = "spearman", adjust = "none", nPerms = 0,
          sigThresh = 1, nPairs = "all"
        ) |>
          suppressMessages()
        as_tibble(observed) |>
          transmute(
            feature = as.character(Gene1), outcome = name,
            rho_first = .data[[paste0(levels(group)[1], "_cor")]],
            rho_second = .data[[paste0(levels(group)[2], "_cor")]],
            z_diff = if_else(is.finite(zScoreDiff), zScoreDiff, NA_real_),
            p = pValDiff, class = Classes
          )
      }) |>
      list_rbind() |>
      mutate(window = window)
  }) |>
    list_rbind()

  tidy <- \(d) {
    d |>
      mutate(
        level = level_name, label = feature_level$label[feature],
        collection = if (is.null(feature_level$collection)) {
          NA_character_
        } else {
          feature_level$collection[feature]
        },
        .before = feature
      ) |>
      mutate(
        fdr = p.adjust(p, "BH"), qvalue = storey_q(p),
        .by = c(level, window, outcome, any_of("group"), collection)
      )
  }
  list(
    within_group = tidy(filter(within_group, !is.na(rho))),
    group_difference = tidy(group_difference),
    flips = tidy(flips)
  )
})
results <- map(set_names(names(per_level[[1]])), \(sheet) list_rbind(map(per_level, sheet)))
names(results$flips)[names(results$flips) == "rho_first"] <- paste0("rho_", group_levels[2])
names(results$flips)[names(results$flips) == "rho_second"] <- paste0("rho_", group_levels[1])

counts <- imap(results, \(d, test) {
  summarise(d,
    n_tested = n(), n_nominal = sum(p < 0.05), n_fdr_05 = sum(fdr < 0.05),
    .by = c(level, window, outcome, any_of("group"))
  ) |>
    mutate(test = test, .before = 1)
}) |>
  list_rbind()

count_plot <- counts |>
  mutate(
    test = if_else(is.na(group), test, paste(test, group)),
    window = factor(window, names(window_of(feature_levels$modules$values)))
  ) |>
  ggplot(aes(n_nominal / n_tested, outcome, fill = test)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_vline(xintercept = 0.05, linetype = "dashed", linewidth = 0.3) +
  facet_grid(level ~ window, scales = "free_y", space = "free_y") +
  scale_x_continuous(labels = scales::label_percent()) +
  labs(
    title = "Share of features at p < 0.05, per test, window and outcome",
    subtitle = sprintf(
      "Groups: %s. Dashed: 5%%. BH within level, window, outcome and collection.",
      paste(group_levels, collapse = " and ")
    ),
    x = "Features at p < 0.05", y = NULL, fill = NULL
  ) +
  theme_bw(base_size = 8) +
  theme(legend.position = "bottom")
pdf(file.path(figures, "06_associate_figures.pdf"), width = 14, height = 10)
print(count_plot)
invisible(dev.off())

write_workbook(
  c(list(counts = counts), map(results, \(d) arrange(d, level, window, outcome, p))),
  c(
    "Tested, p < 0.05 and BH < 0.05 per test, level, window, outcome and group.",
    "Spearman rho and exact p within each group, every feature and outcome.",
    sprintf(
      "limma, feature ~ group * outcome: %s slope minus %s slope, t, p, BH and Storey q.",
      group_levels[2], group_levels[1]
    ),
    "DGCA: each group's Spearman rho, Fisher z difference, its DGCA p, BH, Storey q."
  ),
  file.path(out, "06_associate.xlsx"), inputs, paths
)

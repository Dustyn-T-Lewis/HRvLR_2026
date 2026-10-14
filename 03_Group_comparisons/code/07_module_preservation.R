#' ---
#' title: "Module preservation between HR and LR"
#' params:
#'   grouping: HRvLR
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
#' Do HR and LR share the same co-expression structure? Modules are built inside each group with the
#' settings of `06_Module_Building`, then WGCNA's `modulePreservation()` tests each group's modules
#' in the other group's data (Langfelder et al. 2011). HR/LR only.
#'
#' Two statistics come back per module. Zsummary combines how dense and how connected the module is
#' in the other group, against permutations: above 10 is strong preservation, 2 to 10 weak to
#' moderate, below 2 none (Langfelder et al. 2011). The median rank orders modules against each
#' other and does not depend on module size; lower is better preserved.
#'
#' Each group has 7 or 8 people, so each network rests on 20 to 24 biopsies. The test has little
#' power to show a module is missing: a Zsummary above 2 says the module is found again, not that
#' the two networks are the same.
#'
#' ## Packages

#| label: packages
#| message: false
library(here)
library(dplyr)
library(tibble)
library(tidyr)
library(purrr)
library(ggplot2)
library(WGCNA)
source(here("00_Setup", "Shared_functions", "analysis_reporting_helpers.R"))
disableWGCNAThreads()

#' ## Inputs

#| label: inputs
inputs <- c(
  corrected = "01_Processing/03_Quantification/c_data/DAList_corrected.rds",
  modules = "01_Processing/06_Module_Building/c_data/modules.rds"
)
paths <- map_chr(inputs, here)
stopifnot(file.exists(paths))
corrected <- readRDS(paths[["corrected"]])
modules <- readRDS(paths[["modules"]])
abund <- as.matrix(corrected$data)
meta <- modules$meta
stopifnot(identical(colnames(abund), meta$sample_id))
bicor_args <- list(corType = "bicor", maxPOutliers = 0.05, pearsonFallback = "individual")
permutations <- 200L

#' ## One network per group
#'
#' As in `06_Module_Building`, each group's abundance is centred within person, and people with one
#' biopsy drop out. `apply()` over proteins returns biopsies in rows, the orientation WGCNA expects.

#| label: centre
centred <- function(subjects) {
  kept <- meta |>
    filter(subject %in% subjects) |>
    add_count(subject) |>
    filter(n > 1)
  apply(abund[, kept$sample_id], 1, \(x) x - ave(x, kept$subject))
}
expr <- map(set_names(c("HR", "LR")), \(arm) centred(unique(meta$subject[meta$arm == arm])))
map_int(expr, nrow)

#' `build_arm()` repeats the module construction of `06_Module_Building` on one group's matrix: soft
#' power from `pickSoftThreshold()`, the WGCNA FAQ's signed table as the fallback, and the same
#' `blockwiseModules()` settings and seed.

#| label: build
build_arm <- function(data) {
  sft <- pickSoftThreshold(data,
    powerVector = 1:20, networkType = "signed",
    corFnc = "bicor", corOptions = list(maxPOutliers = 0.05), verbose = 0
  )
  power <- coalesce(
    sft$powerEstimate,
    c(18L, 16L, 14L, 12L)[findInterval(nrow(data), c(20, 31, 41)) + 1]
  )
  net <- do.call(blockwiseModules, c(
    list(data,
      power = power, networkType = "signed", TOMType = "signed",
      deepSplit = 2, minModuleSize = 30, mergeCutHeight = 0.15,
      pamStage = TRUE, pamRespectsDendro = TRUE, maxBlockSize = 2500,
      numericLabels = FALSE, randomSeed = 54321, verbose = 0
    ),
    bicor_args
  ))
  list(power = power, colours = set_names(net$colors, colnames(data)))
}
arm_modules <- map(expr, build_arm)
map(arm_modules, \(x) c(power = x$power, modules = n_distinct(setdiff(x$colours, "grey"))))

#' ## Preservation
#'
#' Each group is the reference once and the test once. `grey` (unassigned) and `gold` (a random
#' sample WGCNA adds as a baseline) are dropped from the table.

#| label: preservation
preservation <- modulePreservation(
  multiData = map(expr, \(data) list(data = data)),
  multiColor = map(arm_modules, "colours"),
  referenceNetworks = 1:2, networkType = "signed",
  corFnc = "bicor", corOptions = "maxPOutliers = 0.05, pearsonFallback = 'individual'",
  nPermutations = permutations, randomSeed = 42, savePermutedStatistics = FALSE,
  verbose = 0
)

#' Each group module is labelled with its best-overlapping module from the full cohort, so the two
#' sets of colour names can be read together.

#| label: table
full <- set_names(modules$membership$module, modules$membership$uniprot_id)
best_match <- imap(arm_modules, \(built, arm) {
  tibble(
    uniprot_id = names(built$colours), arm_module = built$colours, full_module = full[uniprot_id]
  ) |>
    filter(arm_module != "grey") |>
    count(arm_module, full_module) |>
    mutate(share = n / sum(n), .by = arm_module) |>
    slice_max(n, n = 1, by = arm_module, with_ties = FALSE) |>
    transmute(reference = arm, module = arm_module, full_module, full_share = round(share, 2))
}) |>
  list_rbind()
preservation_table <- map(c(HR = "LR", LR = "HR"), \(test) {
  reference <- setdiff(c("HR", "LR"), test)
  pair <- c(paste0("ref.", reference), paste0("inColumnsAlsoPresentIn.", test))
  z <- preservation$preservation$Z[[pair]]
  rank <- preservation$preservation$observed[[pair]]
  tibble(
    reference = reference, test = test, module = rownames(z),
    n_proteins = z$moduleSize, z_summary = z$Zsummary.pres, median_rank = rank$medianRank.pres
  )
}) |>
  list_rbind() |>
  filter(!module %in% c("grey", "gold")) |>
  left_join(best_match, by = c("reference", "module")) |>
  mutate(preservation = cut(
    z_summary, c(-Inf, 2, 10, Inf),
    labels = c("none", "weak to moderate", "strong")
  ))
count(preservation_table, reference, test, preservation)


#| label: fig-preservation
#| fig-cap: "Each group's modules tested in the other; label: best-matching full-cohort module."
#| fig-width: 11
#| fig-height: 7
preservation_table |>
  mutate(direction = paste(reference, "modules tested in", test)) |>
  select(direction, module, n_proteins, full_module,
    `Zsummary` = z_summary, `median rank` = median_rank
  ) |>
  pivot_longer(c(Zsummary, `median rank`), names_to = "statistic") |>
  ggplot(aes(n_proteins, value)) +
  geom_hline(
    data = tibble(statistic = "Zsummary", y = c(2, 10)), aes(yintercept = y),
    linetype = "dashed", colour = "grey55"
  ) +
  geom_point(aes(fill = module), shape = 21, size = 3, colour = "grey30") +
  ggrepel::geom_text_repel(
    aes(label = paste0(module, "\n(", full_module, ")")),
    size = 2.2, lineheight = 0.85, seed = 1, colour = "grey25", min.segment.length = 0
  ) +
  facet_grid(statistic ~ direction, scales = "free_y") +
  scale_fill_identity() +
  scale_x_log10() +
  labs(x = "module size (log scale)", y = NULL) +
  theme_minimal(base_size = 10)

#' ## Validation

#| label: validate
stopifnot(
  setequal(unique(preservation_table$reference), c("HR", "LR")),
  !anyNA(preservation_table$z_summary),
  !anyNA(preservation_table$full_module)
)

#' ## Summary

#| label: summary
bands <- count(preservation_table, reference, preservation)
band <- \(ref, level) sum(filter(bands, reference == ref, preservation == level)$n)

#' HR modules tested in LR: `r band("HR", "strong")` strongly preserved,
#' `r band("HR", "weak to moderate")` weak to moderate, `r band("HR", "none")` not preserved. LR
#' modules tested in HR: `r band("LR", "strong")` strongly preserved,
#' `r band("LR", "weak to moderate")` weak to moderate, `r band("LR", "none")` not preserved.
#'
#' ## Save

#| label: save
out <- here("03_Group_comparisons", params$grouping, "c_data")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write_workbook(
  list(
    preservation = preservation_table,
    arm_networks = tibble(
      arm = names(arm_modules), n_samples = map_int(expr, nrow),
      power = map_int(arm_modules, "power"),
      n_modules = map_int(arm_modules, \(x) n_distinct(setdiff(x$colours, "grey")))
    ),
    arm_membership = imap(arm_modules, \(x, arm) {
      tibble(arm = arm, uniprot_id = names(x$colours), module = x$colours)
    }) |>
      list_rbind()
  ),
  c(
    "Zsummary and median rank of each arm's modules tested in the other arm.",
    "Samples, soft power and module count of each arm's network.",
    "Every protein's module in each arm's network."
  ),
  file.path(out, "07_module_preservation.xlsx"), inputs, paths
)

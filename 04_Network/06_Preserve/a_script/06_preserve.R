# Do HR and LR share the same co-expression structure? Modules are built inside each arm with the
# settings of 01_Build_Modules, then WGCNA's modulePreservation tests each arm's modules in the
# other arm. Zsummary above 10 is strong preservation, 2 to 10 moderate, below 2 none (Langfelder
# et al. 2011); medianRank ranks modules against each other and does not depend on module size.
#
# Each arm has 7 or 8 subjects, so these networks rest on 20 to 24 samples. Read the result as a
# comparison of structure, not as module discovery. The full-cohort modules from 01 are mapped to
# each arm module by overlap so the two sets of names can be read together.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(WGCNA)
})

stage <- here("04_Network", "06_Preserve")
out <- file.path(stage, "c_data")
figure_dir <- file.path(stage, "b_reports")
for (path in c(out, figure_dir)) dir.create(path, recursive = TRUE, showWarnings = FALSE)

inputs <- c(
  imputed = "01_Preprocess/03_Imputation/c_data/DAList_imputed.rds",
  modules = "04_Network/01_Build_Modules/c_data/modules.rds"
)
paths <- map_chr(inputs, here)
if (!all(file.exists(paths))) {
  stop("Run 01_Build_Modules first. Missing: ", paste(inputs[!file.exists(paths)], collapse = ", "))
}
imputed <- readRDS(paths[["imputed"]])
modules <- readRDS(paths[["modules"]])
disableWGCNAThreads()

abund <- as.matrix(imputed$data)
meta <- modules$meta
stopifnot(identical(colnames(abund), meta$sample_id))
bicor_args <- list(corType = "bicor", maxPOutliers = 0.05, pearsonFallback = "individual")
permutations <- 200L

# Subject-centred abundance for one arm, subjects with one biopsy left out as in 01.
centred_arm <- function(arm) {
  kept <- meta |>
    filter(arm == .env$arm) |>
    add_count(subject) |>
    filter(n > 1)
  # apply() over proteins returns samples in rows, the orientation WGCNA expects.
  apply(abund[, kept$sample_id], 1, \(x) x - ave(x, kept$subject))
}
expr <- map(set_names(c("HR", "LR")), centred_arm)
message("samples per arm: ", paste(names(expr), map_int(expr, nrow), collapse = ", "))

# The same construction as 01, arm by arm: pickSoftThreshold's estimate, with the WGCNA FAQ's
# signed table as the fallback power.
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
iwalk(arm_modules, \(x, arm) {
  message(arm, ": power ", x$power, ", ", n_distinct(setdiff(x$colours, "grey")), " modules")
})

preservation <- modulePreservation(
  multiData = map(expr, \(data) list(data = data)),
  multiColor = map(arm_modules, "colours"),
  referenceNetworks = 1:2, networkType = "signed",
  corFnc = "bicor", corOptions = "maxPOutliers = 0.05, pearsonFallback = 'individual'",
  nPermutations = permutations, randomSeed = 42, savePermutedStatistics = FALSE,
  verbose = 0
)

# Each arm module's best-overlapping full-cohort module, for reading the names side by side.
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
    labels = c("none", "moderate", "strong")
  ))
print(as.data.frame(preservation_table), digits = 3)


# ---- figures -------------------------------------------------------------------------------

plot_data <- preservation_table |>
  mutate(direction = paste(reference, "modules tested in", test)) |>
  select(
    direction, module, n_proteins, full_module,
    `Zsummary` = z_summary, `median rank` = median_rank
  ) |>
  pivot_longer(c(Zsummary, `median rank`), names_to = "statistic")
preservation_figure <- ggplot(plot_data, aes(n_proteins, value)) +
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
  labs(
    x = "module size (log scale)", y = NULL, title = "Module preservation between arms",
    subtitle = sprintf(
      "WGCNA modulePreservation, signed bicor, subject-centred, %d permutations", permutations
    ),
    caption = stringr::str_wrap(width = 160, paste(
      "Each point is a module built in one arm and tested in the other; the label gives its",
      "best-overlapping full-cohort module from 01_Build_Modules. Dashed lines mark Zsummary 2",
      "and 10. A lower median rank means stronger preservation. Table:",
      "c_data/06_preserve.xlsx, preservation."
    ))
  ) +
  theme_minimal(base_size = 10) +
  theme(plot.caption = element_text(size = 7, colour = "grey45", hjust = 0))
pdf(file.path(figure_dir, "06_preserve_figures.pdf"), width = 11, height = 8.5)
print(preservation_figure)
invisible(dev.off())

sheets <- list(
  preservation = preservation_table,
  arm_networks = tibble(
    arm = names(arm_modules), n_samples = map_int(expr, nrow),
    power = map_int(arm_modules, "power"),
    n_modules = map_int(arm_modules, \(x) n_distinct(setdiff(x$colours, "grey")))
  ),
  arm_membership = imap(arm_modules, \(x, arm) {
    tibble(arm = arm, uniprot_id = names(x$colours), module = x$colours)
  }) |>
    list_rbind(),
  input_manifest = tibble(
    input = names(inputs), path = unname(inputs), md5 = unname(tools::md5sum(paths))
  ),
  package_versions = sessioninfo::package_info("loaded", dependencies = FALSE) |>
    as_tibble() |>
    select(package, version = loadedversion, source)
)
read_me <- tibble(sheet = names(sheets), holds = c(
  "Zsummary and median rank of each arm's modules tested in the other arm.",
  "Samples, soft power and module count of each arm's network.",
  "Every protein's module in each arm's network.",
  "Files read, with md5.",
  "Packages loaded at run time."
))
writexl::write_xlsx(c(list(read_me = read_me), sheets), file.path(out, "06_preserve.xlsx"))
message("wrote 06_preserve.xlsx and 1 figure")

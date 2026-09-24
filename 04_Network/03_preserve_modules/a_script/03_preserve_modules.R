# Do HR and LR share the same co-expression structure? Modules are built inside each arm with the
# settings of 01_build_modules, then WGCNA's modulePreservation tests each arm's modules in the
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

stage <- here("04_Network", "03_preserve_modules")
out <- file.path(stage, "c_data")
figure_dir <- file.path(stage, "b_reports")
for (path in c(out, figure_dir)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
unlink(list.files(figure_dir, "[.](png|pdf)$", full.names = TRUE))

inputs <- c(
  imputed = "01_Preprocess/03_Imputation/c_data/DAList_imputed.rds",
  modules = "04_Network/01_build_modules/c_data/modules.rds"
)
paths <- map_chr(inputs, here)
if (!all(file.exists(paths))) {
  stop("Run 01_build_modules first. Missing: ", paste(inputs[!file.exists(paths)], collapse = ", "))
}
imputed <- readRDS(paths[["imputed"]])
modules <- readRDS(paths[["modules"]])
manifest <- tibble(
  input = names(inputs), path = unname(inputs), md5 = unname(tools::md5sum(paths))
)
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
map_int(expr, nrow)

# The same construction as 01, arm by arm. The WGCNA FAQ's signed table gives the fallback power.
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
set.seed(42)
arm_modules <- map(expr, build_arm)
map(arm_modules, \(x) c(power = x$power, modules = n_distinct(setdiff(x$colours, "grey"))))

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
  tibble(protein = names(built$colours), arm_module = built$colours, full_module = full[protein]) |>
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
    size = z$moduleSize, z_summary = z$Zsummary.pres, median_rank = rank$medianRank.pres
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

save_figure <- function(figure, name, width, height) {
  walk(c("png", "pdf"), \(extension) {
    ggsave(file.path(figure_dir, paste0(name, ".", extension)), figure,
      width = width, height = height, dpi = 200, bg = "white"
    )
  })
}
plot_data <- preservation_table |>
  mutate(direction = paste(reference, "modules tested in", test)) |>
  select(
    direction, module, size, full_module,
    `Zsummary` = z_summary, `median rank` = median_rank
  ) |>
  pivot_longer(c(Zsummary, `median rank`), names_to = "statistic")
save_figure(
  ggplot(plot_data, aes(size, value)) +
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
        "best-overlapping full-cohort module from 01_build_modules. Dashed lines mark Zsummary 2",
        "and 10. A lower median rank means stronger preservation. Table:",
        "c_data/03_preserve_modules.xlsx, preservation."
      ))
    ) +
    theme_minimal(base_size = 10) +
    theme(plot.caption = element_text(size = 7, colour = "grey45", hjust = 0)),
  "01_preservation", 11, 8
)

packages <- c("here", "WGCNA", "dplyr", "purrr", "ggplot2", "ggrepel")
versions <- tibble(
  package = packages, version = map_chr(packages, \(p) as.character(packageVersion(p)))
)
saveRDS(
  list(
    preservation = preservation_table, arm_modules = map(arm_modules, "colours"),
    powers = map_int(arm_modules, "power"),
    provenance = list(
      created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
      inputs = manifest, packages = versions
    )
  ),
  file.path(out, "module_preservation.rds"),
  compress = "xz"
)
writexl::write_xlsx(
  list(
    preservation = preservation_table,
    arm_membership = imap(arm_modules, \(x, arm) {
      tibble(arm = arm, uniprot_id = names(x$colours), module = x$colours)
    }) |>
      list_rbind(),
    input_manifest = manifest, package_versions = versions
  ),
  file.path(out, "03_preserve_modules.xlsx")
)
combined <- file.path(figure_dir, "03_preserve_modules_figures.pdf")
pages <- setdiff(list.files(figure_dir, "[.]pdf$", full.names = TRUE), combined)
invisible(qpdf::pdf_combine(sort(pages), combined))
message(
  "wrote module_preservation.rds, 03_preserve_modules.xlsx and ", length(pages),
  " figures bundled into ", qpdf::pdf_length(combined), " pages"
)

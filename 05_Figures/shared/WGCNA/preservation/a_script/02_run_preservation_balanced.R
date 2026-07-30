#!/usr/bin/env Rscript
# Balanced variant of 01_run_preservation.R. The unbalanced HR-vs-LR run
# (21 vs 24 samples) found LR preserves into HR slightly better than the
# reverse, but LR's extra samples give it a more precise reference network,
# so that asymmetry could be sample-size artifact rather than biology. This
# script equalizes both arms at 6 complete subjects (18 samples) and repeats
# the LR subject draw 5 times, since which 6 of LR's 8 complete subjects are
# kept is arbitrary.

pacman::p_load(here, dplyr, tibble, purrr, openxlsx, WGCNA, ggplot2, patchwork)
source(here("functions", "shared_wgcna.R"))
source(here("functions", "sweep_drivers.R"))
source(here("functions", "shared_style.R"))

set.seed(42)

imputed <- readRDS(here(
  "02_Normalization", "imputation", "c_data", "DAList_imputed_missforest.rds"
))
stopifnot(identical(imputed$metadata$Col_ID, colnames(imputed$data)))

abund <- as.matrix(imputed$data)
rownames(abund) <- imputed$annotation$gene

meta <- imputed$metadata |>
  transmute(
    sample_id = Col_ID, subject = Subject_ID, group = Group,
    timepoint = Timepoint
  )

mm <- read.xlsx(
  here("05_Figures", "shared", "WGCNA", "c_data", "WGCNA_source_data.xlsx"),
  "module_membership"
)
modules <- setNames(mm$module, mm$gene)[rownames(abund)]
stopifnot(!anyNA(modules))

centred <- centre_within_subject(abund, meta$subject)

subject_timepoints <- meta |>
  group_by(.data$subject, .data$group) |>
  summarise(n_timepoints = dplyr::n_distinct(.data$timepoint), .groups = "drop")

hr_complete <- subject_timepoints |>
  filter(.data$group == "HR", .data$n_timepoints == 3) |>
  pull(.data$subject)
stopifnot(length(hr_complete) == 6L)

lr_complete <- subject_timepoints |>
  filter(.data$group == "LR", .data$n_timepoints == 3) |>
  pull(.data$subject)
stopifnot(length(lr_complete) == 8L)

hr_sample_id <- meta$sample_id[meta$subject %in% hr_complete]
stopifnot(length(hr_sample_id) == 18L)

n_draws <- 5L
lr_draws <- lapply(
  seq_len(n_draws), function(draw) sort(sample(lr_complete, 6))
)

WGCNA::disableWGCNAThreads()

# gold is WGCNA's synthetic random-module benchmark, grey is the unassigned
# background; neither is a reportable co-expression module.
extract_direction <- function(mp, ref_arm, test_arm, draw) {
  z <- mp$preservation$Z[[paste0("ref.", ref_arm)]][[
    paste0("inColumnsAlsoPresentIn.", test_arm)
  ]]
  obs <- mp$preservation$observed[[paste0("ref.", ref_arm)]][[
    paste0("inColumnsAlsoPresentIn.", test_arm)
  ]]
  real_modules <- setdiff(rownames(z), c("grey", "gold"))
  tibble::tibble(
    draw = draw, reference = ref_arm, test = test_arm, module = real_modules,
    n_proteins = obs[real_modules, "moduleSize"],
    Zsummary = z[real_modules, "Zsummary.pres"],
    medianRank = obs[real_modules, "medianRank.pres"]
  )
}

run_draw <- function(draw) {
  lr_subjects <- lr_draws[[draw]]
  lr_sample_id <- meta$sample_id[meta$subject %in% lr_subjects]
  stopifnot(length(lr_sample_id) == 18L)

  mp <- WGCNA::modulePreservation(
    multiData = list(
      HR = list(data = t(centred[, hr_sample_id])),
      LR = list(data = t(centred[, lr_sample_id]))
    ),
    multiColor = list(HR = modules, LR = modules),
    referenceNetworks = c(1, 2),
    nPermutations = 200,
    randomSeed = 42,
    networkType = "signed",
    verbose = 3
  )

  bind_rows(
    extract_direction(mp, "HR", "LR", draw),
    extract_direction(mp, "LR", "HR", draw)
  )
}

preservation <- purrr::map_dfr(seq_len(n_draws), run_draw)

draws <- tibble::tibble(
  draw = rep(seq_len(n_draws), lengths(lr_draws)),
  subject = unlist(lr_draws)
)

design <- purrr::map_dfr(seq_len(n_draws), function(draw) {
  lr_subjects <- lr_draws[[draw]]
  tibble::tibble(
    draw = draw,
    group = c("HR", "LR"),
    n_samples = c(
      length(hr_sample_id), sum(meta$subject %in% lr_subjects)
    ),
    n_subjects = c(length(hr_complete), length(lr_subjects))
  )
})

out_data <- here("05_Figures", "shared", "WGCNA", "preservation", "c_data")
out_reports <- here(
  "05_Figures", "shared", "WGCNA", "preservation", "b_reports"
)
dir.create(out_data, recursive = TRUE, showWarnings = FALSE)
dir.create(out_reports, recursive = TRUE, showWarnings = FALSE)

wb <- createWorkbook()
addWorksheet(wb, "preservation_balanced")
writeData(wb, "preservation_balanced", preservation)
addWorksheet(wb, "draws")
writeData(wb, "draws", draws)
addWorksheet(wb, "design")
writeData(wb, "design", design)
saveWorkbook(
  wb, file.path(out_data, "preservation_balanced.xlsx"),
  overwrite = TRUE
)

threshold <- data.frame(z = c(2, 10), label = c("none", "strong"))

preservation_panel <- function(df, ref_arm, test_arm) {
  d <- df |>
    filter(.data$reference == ref_arm, .data$test == test_arm) |>
    group_by(.data$module) |>
    summarise(
      median_z = stats::median(.data$Zsummary),
      min_z = min(.data$Zsummary),
      max_z = max(.data$Zsummary),
      .groups = "drop"
    )

  ggplot(
    d, aes(stats::reorder(.data$module, -.data$median_z), .data$median_z)
  ) +
    geom_pointrange(
      aes(ymin = .data$min_z, ymax = .data$max_z),
      colour = module_fill(d$module), fatten = 2.5
    ) +
    geom_hline(
      data = threshold, aes(yintercept = .data$z),
      linetype = "dashed", colour = "grey40", inherit.aes = FALSE
    ) +
    geom_text(
      data = threshold, aes(x = Inf, y = .data$z, label = .data$label),
      inherit.aes = FALSE, hjust = 1.05, vjust = -0.4, size = 2.6,
      colour = "grey30"
    ) +
    labs(
      x = NULL, y = "Zsummary.pres (median, range over draws)",
      title = sprintf("reference %s, test %s", ref_arm, test_arm)
    ) +
    FIG_THEME +
    theme(axis.text.x = element_text(angle = 40, hjust = 1))
}

subtitle_lines <- c(
  "Both arms reduced to 6 complete subjects (18 samples each), so",
  "Zsummary is lower than the unbalanced 21-vs-24 run for reasons of",
  "network size alone. Compare the two directions within this figure",
  "only, never against that run."
)

composite <- (
  preservation_panel(preservation, "HR", "LR") |
    preservation_panel(preservation, "LR", "HR")
) +
  plot_annotation(
    title = "Module preservation between HR and LR, balanced subsample",
    subtitle = paste(subtitle_lines, collapse = "\n"),
    theme = theme(
      plot.title = element_text(face = "bold", size = FIG_TITLE_SIZE),
      plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE, color = "grey30")
    )
  )

save_panel(
  composite, file.path(out_reports, "preservation_balanced"),
  width = 240, height = 150
)

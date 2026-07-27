#!/usr/bin/env Rscript
# WGCNA::modulePreservation on the within-subject-centred proteome, split by
# arm after centring. The cohort modules were defined on that centred matrix
# (functions/shared_wgcna.R:63), so preservation is tested on the same
# object a raw-abundance test would not be.

pacman::p_load(here, dplyr, tibble, openxlsx, WGCNA, ggplot2, patchwork)
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
  transmute(sample_id = Col_ID, subject = Subject_ID, group = Group)

mm <- read.xlsx(
  here("04_Figures", "shared", "WGCNA", "c_data", "WGCNA_source_data.xlsx"),
  "module_membership"
)
modules <- setNames(mm$module, mm$gene)[rownames(abund)]
stopifnot(!anyNA(modules))

centred <- centre_within_subject(abund, meta$subject)

hr_id <- meta$sample_id[meta$group == "HR"]
lr_id <- meta$sample_id[meta$group == "LR"]

arm_sizes <- meta |>
  group_by(group) |>
  summarise(
    n_samples = dplyr::n(), n_subjects = dplyr::n_distinct(subject),
    .groups = "drop"
  )

WGCNA::disableWGCNAThreads()

mp <- WGCNA::modulePreservation(
  multiData = list(
    HR = list(data = t(centred[, hr_id])),
    LR = list(data = t(centred[, lr_id]))
  ),
  multiColor = list(HR = modules, LR = modules),
  referenceNetworks = c(1, 2),
  nPermutations = 200,
  randomSeed = 42,
  networkType = "signed",
  verbose = 3
)

# gold is WGCNA's synthetic random-module benchmark, grey is the unassigned
# background; neither is a reportable co-expression module.
extract_direction <- function(mp, ref_arm, test_arm) {
  z <- mp$preservation$Z[[paste0("ref.", ref_arm)]][[
    paste0("inColumnsAlsoPresentIn.", test_arm)
  ]]
  obs <- mp$preservation$observed[[paste0("ref.", ref_arm)]][[
    paste0("inColumnsAlsoPresentIn.", test_arm)
  ]]
  real_modules <- setdiff(rownames(z), c("grey", "gold"))
  tibble::tibble(
    reference = ref_arm, test = test_arm, module = real_modules,
    n_proteins = obs[real_modules, "moduleSize"],
    Zsummary = z[real_modules, "Zsummary.pres"],
    medianRank = obs[real_modules, "medianRank.pres"]
  )
}

preservation <- bind_rows(
  extract_direction(mp, "HR", "LR"),
  extract_direction(mp, "LR", "HR")
)

out_data <- here("04_Figures", "shared", "WGCNA", "preservation", "c_data")
out_reports <- here(
  "04_Figures", "shared", "WGCNA", "preservation", "b_reports"
)
dir.create(out_data, recursive = TRUE, showWarnings = FALSE)
dir.create(out_reports, recursive = TRUE, showWarnings = FALSE)

wb <- createWorkbook()
addWorksheet(wb, "preservation")
writeData(wb, "preservation", preservation)
addWorksheet(wb, "arm_sizes")
writeData(wb, "arm_sizes", arm_sizes)
saveWorkbook(wb, file.path(out_data, "preservation.xlsx"), overwrite = TRUE)

threshold <- data.frame(z = c(2, 10), label = c("none", "strong"))

preservation_panel <- function(df, ref_arm, test_arm) {
  d <- filter(df, .data$reference == ref_arm, .data$test == test_arm)
  ggplot(
    d, aes(stats::reorder(.data$module, -.data$Zsummary), .data$Zsummary)
  ) +
    geom_col(fill = module_fill(d$module), colour = "grey25", linewidth = 0.2) +
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
      x = NULL, y = "Zsummary.pres",
      title = sprintf("reference %s, test %s", ref_arm, test_arm)
    ) +
    FIG_THEME +
    theme(axis.text.x = element_text(angle = 40, hjust = 1))
}

hr_row <- arm_sizes[arm_sizes$group == "HR", ]
lr_row <- arm_sizes[arm_sizes$group == "LR", ]

composite <- (
  preservation_panel(preservation, "HR", "LR") |
    preservation_panel(preservation, "LR", "HR")
) +
  plot_annotation(
    title = "Module preservation between HR and LR",
    subtitle = sprintf(
      paste0(
        "HR %d samples from %d subjects, LR %d samples from %d subjects, ",
        "each arm centred within subject before splitting (one degree of ",
        "freedom spent per subject). At this size a low Zsummary cannot be ",
        "cleanly separated from a genuine difference in module ",
        "architecture; Z > 10 = strong, Z < 2 = none (Langfelder & ",
        "Horvath, 2011)."
      ),
      hr_row$n_samples, hr_row$n_subjects, lr_row$n_samples, lr_row$n_subjects
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = FIG_TITLE_SIZE),
      plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE, color = "grey30")
    )
  )

save_panel(
  composite, file.path(out_reports, "preservation"),
  width = 240, height = 135
)

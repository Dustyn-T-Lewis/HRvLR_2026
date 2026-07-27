#!/usr/bin/env Rscript
# In-sample module scores borrow from the subject they are asked to predict,
# because the cohort network saw that subject. Rebuilding the network inside
# each fold removes the borrowing; the gap between the two is the size of the
# module-definition circularity.

pacman::p_load(here, dplyr, tidyr, readr, openxlsx)
source(here("functions", "loso_wgcna_refit.R"))
source(here("functions", "shared_prediction.R"))
source(here("functions", "pred_features.R"))
source(here("functions", "sweep_grid.R"))

set.seed(54321)

imputed <- readRDS(here(
  "02_Normalization", "imputation", "c_data", "DAList_imputed_missforest.rds"
))
abund <- as.matrix(imputed$data)
rownames(abund) <- imputed$annotation$gene

meta <- imputed$metadata |>
  transmute(
    sample_id = Col_ID, subject = Subject_ID,
    group = factor(Group, levels = c("HR", "LR")),
    timepoint = factor(Timepoint, levels = c("T1", "T2", "T3"))
  )

mm <- read.xlsx(
  here("04_Figures", "shared", "WGCNA", "c_data", "WGCNA_source_data.xlsx"),
  "module_membership"
)
modules_full <- setNames(mm$module, mm$gene)
modules_full <- modules_full[rownames(abund)]

# Subjects are held out by identity, not by constructed sample names: S28 has
# T1 only and S29 has T2 and T3 only, so any paired-name assumption drops them.
subjects <- unique(meta$subject)
folds <- lapply(subjects, function(s) {
  message(sprintf("fold %s", s))
  refit_fold(abund, meta$subject, s, modules_full)
})
names(folds) <- subjects

stability <- bind_rows(lapply(folds, `[[`, "jaccard")) |>
  group_by(.data$full) |>
  summarise(
    mean_jaccard = mean(.data$jaccard),
    min_jaccard = min(.data$jaccard),
    n_folds = dplyr::n(),
    n_failed = sum(.data$jaccard < MIN_JACCARD),
    .groups = "drop"
  ) |>
  rename(module = "full")

# Out-of-fold eigengenes: each subject's samples scored by a network that
# never saw them.
# A fold keeps only the modules whose rebuilt version overlapped the cohort
# module above the Jaccard cut, so fold module sets differ. Binding them raw
# would leave NAs that pred_contrast_matrix passes straight into the models.
# Keep the modules every fold recovered: one a fold could not rebuild without
# a given subject cannot honestly predict that subject either.
complete_modules <- Reduce(intersect, lapply(folds, function(f) {
  names(f$me_hold)
}))
dropped <- setdiff(unique(stability$module), complete_modules)
message(sprintf(
  "modules usable in every fold: %d of %d%s",
  length(complete_modules), length(unique(stability$module)),
  if (length(dropped)) {
    paste0(" (dropped: ", paste(dropped, collapse = ", "), ")")
  } else {
    ""
  }
))

me_refit <- bind_rows(lapply(subjects, function(s) {
  held <- meta$sample_id[meta$subject == s]
  scores <- folds[[s]]$me_hold[complete_modules]
  as.data.frame(scores) |>
    mutate(sample_id = held, .before = 1)
}))
stopifnot(!anyNA(me_refit))

leads <- bind_rows(
  read.xlsx(
    here("04_Figures", "F06_prediction", "MANIFEST.xlsx"),
    "manifest"
  ) |> mutate(stage = "F06_prediction"),
  read.xlsx(
    here("04_Figures", "F05_classification", "MANIFEST.xlsx"),
    "manifest"
  ) |> mutate(stage = "F05_classification")
) |>
  filter(.data$level == "modules", .data$is_lead)

message(sprintf("rescoring %d module leads", nrow(leads)))

bundle <- pred_load()
me_matrix <- me_refit |>
  tibble::column_to_rownames("sample_id") |>
  as.matrix() |>
  t()

rescore <- function(lead) {
  x <- pred_contrast_matrix(me_matrix, bundle$meta, lead$config)
  if (lead$stage == "F06_prediction") {
    al <- align_xy(x, pred_outcome(bundle, lead$phenotype))
    r <- sweep_cont_cell(al$x, al$y, lead$model, lead$phenotype,
      b_grid = c(200L)
    )
    s <- r$summary[r$summary$B == 200, ]
    data.frame(metric_loso_refit = s$q2, perm_p_refit = s$perm_p_q2)
  } else {
    al <- align_xy(x, pred_outcome(bundle, "group"))
    r <- sweep_class_cell(al$x, al$y, lead$model, b_grid = c(200L))
    s <- r$summary[r$summary$B == 200, ]
    data.frame(metric_loso_refit = s$estimate, perm_p_refit = s$perm_p)
  }
}

refit_summary <- bind_rows(lapply(seq_len(nrow(leads)), function(i) {
  lead <- leads[i, ]
  cbind(
    lead[, c("stage", "level", "config", "phenotype", "model")],
    metric_insample = lead$metric_b200,
    rescore(lead)
  )
})) |>
  mutate(drop = .data$metric_insample - .data$metric_loso_refit)

out <- here("04_Figures", "shared", "WGCNA", "loso_refit", "c_data")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write_sweep_workbook(
  file.path(out, "loso_refit.xlsx"),
  list(refit_summary = refit_summary, module_stability = stability)
)

message(sprintf("median drop: %.3f", stats::median(refit_summary$drop)))

pacman::p_load(ggplot2, patchwork)
source(here("functions", "shared_style.R"))
source(here("functions", "sweep_drivers.R"))

slope_long <- refit_summary |>
  mutate(
    cell = sprintf(
      "%s | %s | %s", .data$config, .data$phenotype, .data$model
    )
  ) |>
  pivot_longer(
    c("metric_insample", "metric_loso_refit"),
    names_to = "stage_x", values_to = "metric"
  ) |>
  mutate(
    stage_x = factor(
      ifelse(.data$stage_x == "metric_insample", "in-sample", "LOSO + refit"),
      levels = c("in-sample", "LOSO + refit")
    )
  )

p_slope <- ggplot(slope_long, aes(.data$stage_x, .data$metric,
  group = .data$cell
)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey45") +
  geom_line(colour = "grey40", linewidth = 0.6) +
  geom_point(size = 2, colour = "#B2182B") +
  scale_y_continuous(limits = c(-0.6, 1)) +
  facet_wrap(~ stats::reorder(.data$cell, -.data$metric)) +
  labs(
    x = NULL, y = "metric",
    title = "Module leads, cohort network vs refit inside every fold",
    subtitle = sprintf(
      paste(
        "median drop = %.3f; %d of %d stay above zero.",
        "\nCircularity of module definition only: every",
        "cell\nstill runs on the cohort-imputed proteome."
      ),
      stats::median(refit_summary$drop),
      sum(refit_summary$metric_loso_refit > 0), nrow(refit_summary)
    )
  ) +
  FIG_THEME +
  theme(plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE))

p_stab <- ggplot(stability, aes(
  stats::reorder(.data$module, .data$mean_jaccard), .data$mean_jaccard
)) +
  geom_col(
    fill = module_fill(stability$module), colour = "grey25",
    linewidth = 0.2, width = 0.75
  ) +
  geom_errorbar(aes(ymin = .data$min_jaccard, ymax = .data$mean_jaccard),
    width = 0.25, colour = "grey30"
  ) +
  geom_hline(
    yintercept = MIN_JACCARD, linetype = "dashed",
    colour = "#B2182B"
  ) +
  coord_flip() +
  labs(
    x = NULL, y = "Jaccard overlap with the cohort module",
    subtitle = paste(
      "bar = mean across folds, whisker to the worst.\nDashed line is the 0.5",
      "match cut: only turquoise and\nmagenta clear it in all 16 folds."
    )
  ) +
  FIG_THEME +
  theme(plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE))

save_panel(
  (p_slope | p_stab) + plot_layout(widths = c(1.6, 1)),
  here(
    "04_Figures", "shared", "WGCNA", "loso_refit", "b_reports",
    "loso_refit"
  ),
  width = 260, height = 130
)

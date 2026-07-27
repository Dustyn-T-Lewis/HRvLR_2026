#!/usr/bin/env Rscript
# The preservation asymmetry (LR modules transfer into HR; HR modules do not
# transfer into LR) could reflect HR subjects resembling each other less
# than LR subjects do. This measures that directly on the same centred,
# non-grey protein space the modules were defined on: each subject's own
# protein-by-protein correlation matrix, compared pairwise within an arm by
# the Spearman correlation of the vectorised upper triangles.

pacman::p_load(here, dplyr, tibble, openxlsx, ggplot2)
source(here("functions", "shared_wgcna.R"))
source(here("functions", "shared_style.R"))

set.seed(42)
n_permutations <- 2000

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
genes_keep <- names(modules)[modules != "grey"]

centred <- centre_within_subject(abund, meta$subject)

timepoint_counts <- meta |>
  count(subject, group, name = "n_timepoints")

dropped <- filter(timepoint_counts, n_timepoints < 3)
kept <- filter(timepoint_counts, n_timepoints >= 3)

if (nrow(dropped)) {
  detail <- paste(
    sprintf("%s (%d)", dropped$subject, dropped$n_timepoints),
    collapse = ", "
  )
  message(sprintf(
    paste0(
      "dropped %s: fewer than 3 timepoints cannot support a subject-level ",
      "correlation matrix"
    ),
    detail
  ))
}
message(sprintf(
  "kept %d HR and %d LR subjects (>=3 timepoints each)",
  sum(kept$group == "HR"), sum(kept$group == "LR")
))

subjects_kept <- kept$subject
arm_lookup <- setNames(kept$group, kept$subject)

corr_by_subject <- lapply(subjects_kept, function(s) {
  cols <- meta$sample_id[meta$subject == s]
  stats::cor(t(centred[genes_keep, cols]), method = "pearson")
})
names(corr_by_subject) <- subjects_kept

upper_idx <- upper.tri(corr_by_subject[[1]])
upper_vecs <- lapply(corr_by_subject, function(m) m[upper_idx])

subject_pairs <- utils::combn(subjects_kept, 2, simplify = FALSE)

full_pairwise <- bind_rows(lapply(subject_pairs, function(p) {
  tibble(
    subject_a = p[1], subject_b = p[2],
    arm_a = arm_lookup[[p[1]]], arm_b = arm_lookup[[p[2]]],
    similarity = stats::cor(
      upper_vecs[[p[1]]], upper_vecs[[p[2]]],
      method = "spearman", use = "pairwise.complete.obs"
    )
  )
}))

pairwise <- full_pairwise |>
  filter(arm_a == arm_b) |>
  transmute(arm = arm_a, subject_a, subject_b, similarity)

summary_tbl <- pairwise |>
  group_by(arm) |>
  summarise(
    n_subjects = n_distinct(c(subject_a, subject_b)),
    n_pairs = n(),
    median = stats::median(similarity),
    iqr_lo = unname(stats::quantile(similarity, 0.25)),
    iqr_hi = unname(stats::quantile(similarity, 0.75)),
    min = min(similarity),
    max = max(similarity),
    .groups = "drop"
  )

group_median_diff <- function(labels) {
  arm_a <- labels[full_pairwise$subject_a]
  arm_b <- labels[full_pairwise$subject_b]
  same_arm <- arm_a == arm_b
  med <- tapply(full_pairwise$similarity[same_arm], arm_a[same_arm], median)
  unname(med["HR"] - med["LR"])
}

observed_diff <- group_median_diff(arm_lookup)
perm_diffs <- vapply(seq_len(n_permutations), function(i) {
  shuffled <- setNames(sample(unname(arm_lookup)), subjects_kept)
  group_median_diff(shuffled)
}, numeric(1))
perm_p <- (sum(abs(perm_diffs) >= abs(observed_diff)) + 1) /
  (n_permutations + 1)

test_tbl <- tibble(
  observed_median_diff = observed_diff,
  permutation_p = perm_p,
  n_permutations = n_permutations,
  n_dropped_subjects = nrow(dropped),
  dropped_subjects = paste(dropped$subject, collapse = "; "),
  dropped_timepoints = paste(dropped$n_timepoints, collapse = "; ")
)

out_data <- here("04_Figures", "shared", "WGCNA", "heterogeneity", "c_data")
out_reports <- here(
  "04_Figures", "shared", "WGCNA", "heterogeneity", "b_reports"
)
dir.create(out_data, recursive = TRUE, showWarnings = FALSE)
dir.create(out_reports, recursive = TRUE, showWarnings = FALSE)

wb <- createWorkbook()
addWorksheet(wb, "pairwise")
writeData(wb, "pairwise", pairwise)
addWorksheet(wb, "summary")
writeData(wb, "summary", summary_tbl)
addWorksheet(wb, "test")
writeData(wb, "test", test_tbl)
saveWorkbook(wb, file.path(out_data, "heterogeneity.xlsx"), overwrite = TRUE)

hr_row <- summary_tbl[summary_tbl$arm == "HR", ]
lr_row <- summary_tbl[summary_tbl$arm == "LR", ]

heterogeneity_plot <- ggplot(
  pairwise |> mutate(arm = factor(arm, levels = c("HR", "LR"))),
  aes(arm, similarity, colour = arm)
) +
  geom_boxplot(width = 0.35, outlier.shape = NA, colour = "grey30", fill = NA) +
  geom_jitter(width = 0.12, height = 0, size = 1.8, alpha = 0.85) +
  scale_colour_manual(values = GROUP_COLORS, guide = "none") +
  labs(
    x = NULL, y = "pairwise similarity (Spearman rho of upper triangles)",
    title = "Within-arm subject heterogeneity in protein correlation structure",
    subtitle = sprintf(
      paste0(
        "HR %d subjects, %d pairs; LR %d subjects, %d pairs. Pairs share ",
        "subjects and are not independent, so significance comes from a\n",
        "%d-permutation shuffle of arm labels over subjects, not a rank ",
        "test on the pairs. Observed median difference %.3f, p = %.3f."
      ),
      hr_row$n_subjects, hr_row$n_pairs, lr_row$n_subjects, lr_row$n_pairs,
      n_permutations, observed_diff, perm_p
    )
  ) +
  FIG_THEME

save_panel(
  heterogeneity_plot, file.path(out_reports, "heterogeneity"),
  width = 150, height = 115
)

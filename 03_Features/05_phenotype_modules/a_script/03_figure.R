#!/usr/bin/env Rscript
# One composite of the label-free phenotype-association stage: the
# module map, the pathway map, the delta clusters, and the single
# near-miss shown at full size. Reads only this stage's c_data.
#
# Stars mark empirical p < 0.05 within one timepoint; black boxes mark
# the pre-declared consistency call (same sign at all three timepoints,
# p < 0.05 at two or more). The consistency counts' own permutation
# verdicts sit in the panel subtitles — that second-level null is the
# panel's actual result, and nothing here beats it.

pacman::p_load(
  here, dplyr, tidyr, readr, stringr, purrr, ggplot2, patchwork
)
source(here("functions", "shared_style.R"))

DAT <- here("03_Features", "05_phenotype_modules", "c_data")
OUT <- here("03_Features", "05_phenotype_modules", "b_reports")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

PHENO_LABELS <- c(
  comp_hypertrophy = "Composite", d_fcsa_I = "Δ fCSA I",
  d_fcsa_II = "Δ fCSA II", d_mcsa = "Δ mCSA",
  d_1rm_legpress = "Δ 1RM press", d_1rm_ext = "Δ 1RM ext"
)
REPRODUCIBLE <- c("pink", "turquoise")

maps <- bind_rows(
  read_csv(file.path(DAT, "01_module_trait.csv"), show_col_types = FALSE),
  read_csv(file.path(DAT, "01_pathway_trait.csv"), show_col_types = FALSE)
) |>
  mutate(
    pheno_label = factor(PHENO_LABELS[.data$phenotype],
      levels = PHENO_LABELS
    )
  ) |>
  mutate(
    consistent = n_distinct(sign(.data$rho)) == 1 &
      sum(.data$emp_p < 0.05) >= 2,
    .by = c("level", "feature", "phenotype")
  )

consistency <- read_csv(
  file.path(DAT, "01_consistency.csv"),
  show_col_types = FALSE
)

fill_rho <- scale_fill_gradient2(
  low = DIR_COLORS[["Down"]], mid = "white", high = DIR_COLORS[["Up"]],
  limits = c(-1, 1), name = "Spearman ρ"
)

trait_heatmap <- function(df) {
  ggplot(df, aes(.data$pheno_label, .data$feature_label)) +
    geom_tile(aes(fill = .data$rho), colour = "grey85", linewidth = 0.2) +
    geom_tile(
      data = filter(df, .data$consistent),
      fill = NA, colour = "black", linewidth = 0.5
    ) +
    geom_text(
      data = filter(df, .data$emp_p < 0.05), label = "*",
      size = 3, vjust = 0.75
    ) +
    facet_wrap(~timepoint, nrow = 1) +
    fill_rho +
    labs(x = NULL, y = NULL) +
    FIG_THEME +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )
}

modules <- maps |>
  filter(.data$level == "module") |>
  mutate(feature_label = ifelse(
    .data$feature %in% REPRODUCIBLE,
    paste0(.data$feature, " (reproducible)"), .data$feature
  ))

cons_line <- function(lvl) {
  rows <- filter(consistency, .data$level == lvl, .data$n_observed > 0)
  if (!nrow(rows)) {
    return("no phenotype has a consistent feature")
  }
  paste(sprintf(
    "%s %d (p = %.2f)",
    PHENO_LABELS[rows$phenotype], rows$n_observed, rows$emp_p
  ), collapse = ", ")
}

p_a <- trait_heatmap(modules) +
  labs(
    title = "WGCNA eigengenes vs phenotypes, per timepoint",
    subtitle = paste0(
      "consistency counts vs their own null — ", cons_line("module")
    )
  )

pathway_keep <- maps |>
  filter(.data$level == "pathway") |>
  mutate(
    keep = sum(.data$emp_p < 0.05) >= 2,
    .by = c("feature", "phenotype")
  ) |>
  mutate(keep_set = any(.data$keep | .data$consistent), .by = "feature") |>
  filter(.data$keep_set) |>
  mutate(feature_label = clean_pathway_name(.data$feature, 32))

p_b <- trait_heatmap(pathway_keep) +
  labs(
    title = sprintf(
      "Hallmark singscores, p < 0.05 twice for one phenotype (%d of 57)",
      n_distinct(pathway_keep$feature)
    ),
    subtitle = paste0("consistency vs null — ", cons_line("pathway"))
  )

clusters <- read_csv(file.path(DAT, "02_clusters.csv"),
  show_col_types = FALSE
)
cluster_trait <- read_csv(file.path(DAT, "02_cluster_trait.csv"),
  show_col_types = FALSE
) |>
  mutate(
    pheno_label = factor(PHENO_LABELS[.data$phenotype],
      levels = PHENO_LABELS
    ),
    feature_label = sprintf(
      "%s (n = %d)", .data$feature,
      c(table(clusters$cluster))[.data$feature]
    )
  )

p_c <- ggplot(cluster_trait, aes(.data$pheno_label, .data$feature_label)) +
  geom_tile(aes(fill = .data$rho), colour = "grey85", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.2f", .data$rho)), size = 2.6) +
  fill_rho +
  guides(fill = "none") +
  labs(
    x = NULL, y = NULL,
    title = "Training-response delta clusters (k = 2)",
    subtitle = "no cluster reaches permutation\np < 0.05 on any phenotype"
  ) +
  FIG_THEME +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

eig <- read_csv(
  here(
    "03_Analysis", "categorical", "03_WGCNA", "c_data", "wgcna_eigengene.csv"
  ),
  show_col_types = FALSE
) |>
  filter(.data$group_id == "brown") |>
  mutate(
    subject = sub("_T\\d$", "", .data$sample_id),
    timepoint = sub("^.*_", "", .data$sample_id)
  )
pheno <- read_csv(
  here("00_input", "c_data", "phenotype.csv"),
  show_col_types = FALSE
) |>
  mutate(subject = sub("^(HR|LR)_", "", .data$subject))
brown <- inner_join(eig, pheno, by = "subject")
brown_rho <- maps |>
  filter(.data$feature == "brown", .data$phenotype == "d_mcsa") |>
  select("timepoint", "rho", "emp_p")

p_d <- ggplot(brown, aes(.data$d_mcsa, .data$ME)) +
  geom_point(aes(colour = .data$group_arm), size = 1.8) +
  geom_smooth(
    method = "lm", formula = y ~ x, se = FALSE,
    colour = "grey40", linewidth = 0.5
  ) +
  geom_text(
    data = brown_rho,
    aes(label = sprintf("ρ = %.2f, p = %.3f", .data$rho, .data$emp_p)),
    x = Inf, y = -Inf, hjust = 1.1, vjust = -0.8, size = 2.7,
    inherit.aes = FALSE
  ) +
  facet_wrap(~timepoint, nrow = 1) +
  scale_colour_manual(values = GROUP_COLORS, name = NULL) +
  labs(
    x = "Δ mCSA (cm²)", y = "brown eigengene",
    title = "The near-miss: brown vs Δ mCSA",
    subtitle = paste(
      "stable sign, consistency p = 0.057;\nbrown fails",
      "LOSO network rebuilding"
    )
  ) +
  FIG_THEME +
  theme(legend.position = "right")

fig <- p_a / p_b / ((p_c | p_d) + plot_layout(widths = c(1, 1.5))) +
  plot_layout(heights = c(1.15, 0.85, 0.9)) +
  plot_annotation(
    tag_levels = "A",
    title = paste(
      "Label-free feature clusters vs continuous phenotypes,",
      "tested independently per timepoint"
    ),
    subtitle = paste(
      "14 subjects with all three timepoints; B = 1000",
      "subject-permutation nulls; timepoints share subjects, so",
      "consistency is stability, not replication"
    ),
    theme = theme(
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 9, colour = "grey30")
    )
  )

save_panel(
  fig, file.path(OUT, "F_phenotype_modules"),
  width = 260, height = 300
)
message("wrote ", file.path(OUT, "F_phenotype_modules.{pdf,png}"))

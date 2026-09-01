#!/usr/bin/env Rscript
# Composite for the change configurations and the pooled response: what
# training and the acute bout do on average, and whether any change
# tracks the phenotypes at module or pathway level. Reads only this
# stage's c_data.
#
# Panel C carries the engine disagreement on purpose: fgsea NES on the
# x-axis, black boxes only where limma::fry confirms at FDR < 0.05.
# fgsea permutes gene labels and inflates co-regulated sets (OxPhos most
# of all in muscle); fry rotates residuals and keeps the correlation.
# Panel D applies the same scepticism to the phenotype rankings, where
# no fry counterpart exists: the recurring OxPhos rows flip sign
# between the two strength measures inside one config, which is the
# signature of a correlated-gene artifact, not biology.

pacman::p_load(
  here, dplyr, tidyr, readr, stringr, purrr, ggplot2, ggrepel, patchwork
)
source(here("functions", "shared_style.R"))

DAT <- here("03_Features", "supplementary", "phenotype_modules", "c_data")
OUT <- here("03_Features", "supplementary", "phenotype_modules", "b_reports")

PHENO_LABELS <- c(
  comp_hypertrophy = "Composite", d_fcsa_I = "Δ fCSA I",
  d_fcsa_II = "Δ fCSA II", d_mcsa = "Δ mCSA",
  d_1rm_legpress = "Δ 1RM press", d_1rm_ext = "Δ 1RM ext"
)

fill_rho <- scale_fill_gradient2(
  low = DIR_COLORS[["Down"]], mid = "white", high = DIR_COLORS[["Up"]],
  limits = c(-1, 1), name = "Spearman ρ"
)

change <- read_csv(file.path(DAT, "04_change_trait.csv"),
  show_col_types = FALSE
) |>
  mutate(
    pheno_label = factor(PHENO_LABELS[.data$phenotype],
      levels = PHENO_LABELS
    ),
    config = factor(.data$config, levels = c("training", "acute"))
  )

p_a <- change |>
  filter(.data$level == "module") |>
  ggplot(aes(.data$pheno_label, .data$feature)) +
  geom_tile(aes(fill = .data$rho), colour = "grey85", linewidth = 0.2) +
  geom_text(
    data = \(d) filter(d, .data$emp_p < 0.05), label = "*",
    size = 3, vjust = 0.75
  ) +
  facet_wrap(~config, nrow = 1) +
  fill_rho +
  labs(
    x = NULL, y = NULL,
    title = "Module eigengene changes vs phenotypes",
    subtitle = "3 nominal hits of 144 tests;\nabout 7 expected by chance"
  ) +
  FIG_THEME +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pooled <- read_csv(file.path(DAT, "05_pooled_dep.csv"),
  show_col_types = FALSE
) |>
  filter(!is.na(.data$p)) |>
  mutate(contrast = factor(.data$contrast,
    levels = c("Training_All", "Acute_All")
  ))
counts <- pooled |>
  summarise(bh05 = sum(.data$bh < 0.05), .by = "contrast")
top_genes <- pooled |>
  filter(.data$bh < 0.05) |>
  slice_min(.data$bh, n = 6, by = "contrast")

p_b <- ggplot(pooled, aes(.data$logFC, -log10(.data$p))) +
  geom_point(size = 0.4, alpha = 0.25, colour = "grey40") +
  geom_point(
    data = filter(pooled, .data$bh < 0.05),
    size = 0.7, colour = GROUP_COLORS[["HR"]]
  ) +
  geom_text_repel(
    data = top_genes, aes(label = .data$gene),
    size = 2.4, max.overlaps = 20, seed = 42
  ) +
  facet_wrap(~contrast, nrow = 1) +
  labs(
    x = expression(log[2] ~ "fold change"),
    y = expression(-log[10] ~ p),
    title = "Pooled all-subject contrasts (n = 16)",
    subtitle = sprintf(
      "BH < 0.05: %s;\nblue = significant",
      paste(sprintf(
        "%s %d", counts$contrast, counts$bh05
      ), collapse = ", ")
    )
  ) +
  FIG_THEME

pooled_pw <- read_csv(file.path(DAT, "05_pooled_pathways.csv"),
  show_col_types = FALSE
) |>
  filter(!is.na(.data$fgsea_nes)) |>
  mutate(
    fry_sig = !is.na(.data$fry_fdr) & .data$fry_fdr < 0.05,
    pathway_label = clean_pathway_name(.data$pathway, 34),
    contrast = factor(.data$contrast,
      levels = c("Training_All", "Acute_All")
    )
  ) |>
  filter(.data$fry_sig | .data$fgsea_padj < 0.01)

p_c <- ggplot(pooled_pw, aes(
  .data$fgsea_nes, reorder(.data$pathway_label, .data$fgsea_nes)
)) +
  geom_vline(xintercept = 0, colour = "grey70", linetype = 2) +
  geom_segment(aes(x = 0, xend = .data$fgsea_nes), colour = "grey75") +
  geom_point(aes(colour = .data$fgsea_nes > 0, shape = .data$fry_sig),
    size = 2.4
  ) +
  scale_colour_manual(
    values = c(`TRUE` = DIR_COLORS[["Up"]], `FALSE` = DIR_COLORS[["Down"]]),
    guide = "none"
  ) +
  scale_shape_manual(
    values = c(`TRUE` = 16, `FALSE` = 1),
    labels = c(`TRUE` = "fry FDR < 0.05", `FALSE` = "fgsea only"),
    name = NULL
  ) +
  facet_wrap(~contrast, nrow = 1, scales = "free_y") +
  labs(
    x = "fgsea NES", y = NULL,
    title = "Pooled-contrast pathways: two engines",
    subtitle = paste(
      "filled = rotation-confirmed;\nopen = fgsea-only,",
      "incl. OxPhos at padj 4e-17 that fry rejects"
    )
  ) +
  FIG_THEME

fp <- read_csv(file.path(DAT, "05_fgsea_pheno.csv"),
  show_col_types = FALSE
) |>
  mutate(
    pheno_label = factor(PHENO_LABELS[.data$phenotype],
      levels = PHENO_LABELS
    ),
    config = factor(.data$config, levels = c("training", "acute"))
  ) |>
  mutate(keep = any(.data$padj < 0.05), .by = c("pathway", "config")) |>
  filter(.data$keep) |>
  mutate(pathway_label = clean_pathway_name(.data$pathway, 30))

p_d <- ggplot(fp, aes(.data$pheno_label, .data$pathway_label)) +
  geom_tile(aes(fill = .data$NES), colour = "grey85", linewidth = 0.2) +
  geom_text(
    data = filter(fp, .data$padj < 0.05), label = "*",
    size = 3, vjust = 0.75
  ) +
  facet_wrap(~config, nrow = 1, scales = "free_y") +
  scale_fill_gradient2(
    low = DIR_COLORS[["Down"]], mid = "white", high = DIR_COLORS[["Up"]],
    name = "fgsea NES"
  ) +
  labs(
    x = NULL, y = NULL,
    title = "fgsea on phenotype-correlation rankings",
    subtitle = paste(
      "gene-permutation null only;\nOxPhos flips sign between",
      "the strength measures in training"
    )
  ) +
  FIG_THEME +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

fig <- (p_a | p_b) / (p_c | p_d) +
  plot_layout(heights = c(0.95, 1.05)) +
  plot_annotation(
    tag_levels = "A",
    title = paste(
      "Training and acute change configurations,",
      "pooled and per phenotype"
    ),
    subtitle = paste(
      "Per-subject deltas need no design contrast; the pooled",
      "Training_All / Acute_All contrasts supply the average-response",
      "backdrop the nine canonical contrasts never tested"
    ),
    theme = theme(
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 9, colour = "grey30")
    )
  )

save_panel(
  fig, file.path(OUT, "F_change_response"),
  width = 330, height = 240
)
message("wrote ", file.path(OUT, "F_change_response.{pdf,png}"))

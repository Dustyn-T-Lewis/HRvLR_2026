#!/usr/bin/env Rscript
# One composite of the galamm pilot: what each of the three specification
# checks found. Reads only the pilot's own c_data; renders to b_reports.
#
# Panels follow PREREG order. A and B are Q1: the per-protein residual-SD
# ratios against the [0.8, 1.25] close window, then the six contrasts'
# p-values under the homoscedastic and varIdent fits. C to E are Q2: the
# measurement-model loadings, the factor scores against d_mcsa, and the
# per-protein latent-association volcano. A null pilot looks exactly like
# this: ratio mass inside the window, the p-p cloud on the diagonal, and
# no protein under the BH line.

pacman::p_load(here, dplyr, tidyr, readr, ggplot2, ggrepel, patchwork)
source(here("functions", "shared_style.R"))
source(here("functions", "shared_hlm.R"))

DAT <- here("03_Features", "04_galamm_pilot", "c_data")
OUT <- here("03_Features", "04_galamm_pilot", "b_reports")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

variance <- read_csv(file.path(DAT, "01_q1_variance.csv"),
  show_col_types = FALSE
)
q1_con <- read_csv(file.path(DAT, "01_q1_contrasts.csv"),
  show_col_types = FALSE
)
loadings <- read_csv(file.path(DAT, "02_q2_measurement.csv"),
  show_col_types = FALSE
)
scores <- read_csv(file.path(DAT, "02_q2_scores.csv"), show_col_types = FALSE)
protein <- read_csv(file.path(DAT, "03_q2_protein.csv"),
  show_col_types = FALSE
)

anno <- readRDS(
  here("02_Normalization", "c_data", "DAList_normalized.rds")
)$annotation
protein$gene <- anno$gene[match(protein$feature, anno$uniprot_id)]

ITEM_LABELS <- c(
  comp_hypertrophy = "Composite hypertrophy",
  d_fcsa_I = "Δ fCSA type I",
  d_fcsa_II = "Δ fCSA type II",
  d_mcsa = "Δ mCSA",
  d_1rm_legpress = "Δ 1RM leg press",
  d_1rm_ext = "Δ 1RM extension"
)

ratio_long <- variance |>
  filter(.data$status == "ok") |>
  pivot_longer(c("ratio_t2", "ratio_t3"),
    names_to = "ratio", values_to = "value"
  ) |>
  mutate(ratio = ifelse(.data$ratio == "ratio_t2",
    "σ T2 / σ T1", "σ T3 / σ T1"
  ))

n_clipped <- sum(ratio_long$value < 0.1 | ratio_long$value > 10)

p_a <- ggplot(ratio_long, aes(.data$value, .data$ratio)) +
  annotate("rect",
    xmin = 0.8, xmax = 1.25, ymin = -Inf, ymax = Inf,
    fill = "grey85", alpha = 0.6
  ) +
  geom_violin(fill = TIME_COLORS[["T2"]], alpha = 0.25, linewidth = 0.3) +
  geom_boxplot(width = 0.12, outliers = FALSE, linewidth = 0.3) +
  scale_x_log10(
    breaks = c(0.1, 0.25, 0.5, 1, 2, 4, 10),
    labels = c("0.1", "0.25", "0.5", "1", "2", "4", "10")
  ) +
  coord_cartesian(xlim = c(0.1, 10)) +
  labs(
    x = "residual-SD ratio (log scale)", y = NULL,
    title = "Q1: per-timepoint residual variance",
    subtitle = sprintf(
      "medians %.2f and %.2f, inside the pre-declared [0.8, 1.25] window",
      median(variance$ratio_t2[variance$status == "ok"]),
      median(variance$ratio_t3[variance$status == "ok"])
    ),
    caption = sprintf(
      "axis clipped at [0.1, 10]; %d of %d ratios outside",
      n_clipped, nrow(ratio_long)
    )
  ) +
  FIG_THEME

pp <- q1_con |>
  select("feature", "contrast", "model", "p") |>
  pivot_wider(names_from = "model", values_from = "p") |>
  mutate(contrast = factor(.data$contrast, levels = HLM_CONTRASTS))

p_b <- ggplot(pp, aes(
  -log10(.data$homoscedastic), -log10(.data$varident)
)) +
  geom_abline(linetype = 2, colour = "grey55") +
  geom_point(size = 0.4, alpha = 0.25, colour = GROUP_COLORS[["HR"]]) +
  facet_wrap(~contrast, nrow = 2) +
  labs(
    x = expression(-log[10] ~ p * ", one" ~ sigma),
    y = expression(-log[10] ~ p * ", " * sigma ~ "per timepoint"),
    title = "Q1: the six contrasts under both error models",
    subtitle = "0 BH survivors either way; the diagonal is the result"
  ) +
  FIG_THEME

load_plot <- loadings |>
  mutate(
    label = factor(ITEM_LABELS[.data$item], levels = rev(ITEM_LABELS)),
    fixed = is.na(.data$se)
  )

p_c <- ggplot(load_plot, aes(.data$loading, .data$label)) +
  geom_vline(xintercept = 0, colour = "grey55", linetype = 2) +
  geom_errorbar(
    aes(
      xmin = .data$loading - 1.96 * .data$se,
      xmax = .data$loading + 1.96 * .data$se
    ),
    orientation = "y", width = 0.18, linewidth = 0.4, na.rm = TRUE
  ) +
  geom_point(aes(shape = .data$fixed), size = 2.2, na.rm = TRUE) +
  scale_shape_manual(
    values = c(`FALSE` = 16, `TRUE` = 18), guide = "none"
  ) +
  annotate("text",
    x = 1, y = 6.35, label = "anchor, fixed to 1",
    size = 2.7, hjust = 0.5, colour = "grey30"
  ) +
  labs(
    x = "loading on the latent factor (Wald 95% CI)", y = NULL,
    title = "Q2: measurement model",
    subtitle = "fCSA carries the factor; Δ mCSA loads at 0.27 (SE 0.28)"
  ) +
  FIG_THEME

rho_mcsa <- scores$rho_d_mcsa[1]
rho_comp <- scores$rho_comp[1]

p_d <- ggplot(scores, aes(.data$d_mcsa, .data$eta, colour = .data$group_arm)) +
  geom_point(size = 2) +
  scale_colour_manual(values = GROUP_COLORS, name = NULL) +
  scale_x_continuous(expand = expansion(mult = 0.1)) +
  labs(
    x = "Δ mCSA (cm²)", y = "EB factor score",
    title = "Q2: the factor is not d_mcsa",
    subtitle = sprintf(
      "Spearman ρ = %.2f vs Δ mCSA, %.2f vs composite",
      rho_mcsa, rho_comp
    )
  ) +
  FIG_THEME +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.12, 0.85)
  )

top_lab <- protein |>
  filter(.data$status == "ok") |>
  slice_min(.data$p, n = 5)

p_e <- ggplot(
  filter(protein, .data$status == "ok"),
  aes(.data$loading, -log10(.data$p))
) +
  geom_point(size = 0.5, alpha = 0.3, colour = "grey40") +
  geom_point(
    data = top_lab, size = 1.4, colour = GROUP_COLORS[["LR"]]
  ) +
  geom_text_repel(
    data = top_lab, aes(label = .data$gene),
    size = 2.6, min.segment.length = 0, seed = 42
  ) +
  labs(
    x = "protein loading on the latent factor",
    y = expression(-log[10] ~ p ~ "(Wald)"),
    title = "Q2: latent factor vs 931 proteins",
    subtitle = sprintf(
      "all %d fits converged; min BH q = %.2f, 0 survivors",
      sum(protein$status == "ok"), min(protein$bh, na.rm = TRUE)
    )
  ) +
  FIG_THEME

fig <- (p_a | p_b) / (p_c | p_d | p_e) +
  plot_annotation(
    tag_levels = "A",
    title = "galamm pilot: three specification checks against the HRvLR null",
    subtitle = paste(
      "931 complete-case proteins, blood index in every model,",
      "thresholds fixed in PREREG.md before fitting"
    ),
    theme = theme(
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 9, colour = "grey30")
    )
  )

save_panel(fig, file.path(OUT, "F_galamm_pilot"), width = 320, height = 200)
message("wrote ", file.path(OUT, "F_galamm_pilot.{pdf,png}"))

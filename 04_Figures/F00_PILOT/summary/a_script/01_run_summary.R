#!/usr/bin/env Rscript
# F00_PILOT summary: how much the proteome moves (magnitude) against how much HR and
# LR move together (concordance), for training and the acute bout side by side. The
# five-panel composites carry the per-condition detail; this is the one read-out that
# puts magnitude and concordance in the same frame before either figure is finalised.
pacman::p_load(here, dplyr, tidyr, readr, ggplot2, patchwork)

source(here("04_Figures", "functions", "shared_style.R"))

unit <- here("04_Figures", "F00_PILOT", "summary")
rpt <- file.path(unit, "b_reports")
dat <- file.path(unit, "c_data")

dep <- read_csv(
  here("03_DEP", "a_non_imputed", "c_data", "03_combined_results.csv"),
  show_col_types = FALSE
)

conditions <- tibble::tribble(
  ~phase,      ~c_hi,         ~c_lo,         ~c_int,
  "Training",  "Training_HR", "Training_LR", "Training_Interaction",
  "Acute",     "Acute_HR",    "Acute_LR",    "Acute_Interaction"
)

per_protein <- conditions |>
  rowwise() |>
  reframe(
    phase = phase,
    gene = dep$gene,
    lfc_HR = dep[[paste0("logFC_", c_hi)]],
    lfc_LR = dep[[paste0("logFC_", c_lo)]],
    sig_int = dep[[paste0("sig_pi_", c_int)]] != 0
  ) |>
  filter(!is.na(lfc_HR), !is.na(lfc_LR)) |>
  mutate(phase = factor(phase, levels = c("Training", "Acute")))

stopifnot(nrow(per_protein) > 0, all(c("Training", "Acute") %in% per_protein$phase))

concordance <- per_protein |>
  group_by(phase) |>
  summarise(
    n = n(),
    rho = cor(lfc_HR, lfc_LR, method = "spearman"),
    same_dir = mean(sign(lfc_HR) == sign(lfc_LR)),
    n_interaction = sum(sig_int),
    .groups = "drop"
  ) |>
  mutate(
    ci = lapply(seq_len(n()), function(i) fisher_z_ci(rho[i], n[i])),
    rho_lo = vapply(ci, function(z) z[["lo"]], numeric(1)),
    rho_hi = vapply(ci, function(z) z[["hi"]], numeric(1))
  ) |>
  select(-ci)

magnitude <- per_protein |>
  pivot_longer(c(lfc_HR, lfc_LR), names_to = "group", values_to = "lfc") |>
  mutate(group = sub("^lfc_", "", group), abs_lfc = abs(lfc)) |>
  group_by(phase, group) |>
  summarise(
    med_abs = median(abs_lfc),
    p90_abs = quantile(abs_lfc, 0.90),
    .groups = "drop"
  )

# Panel A: median absolute logFC per group, training vs acute. The proteome moves a
# little more in HR than LR in both conditions, and a little more in the acute bout.
pa <- ggplot(magnitude, aes(phase, med_abs, fill = group)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_errorbar(
    aes(ymin = med_abs, ymax = p90_abs),
    position = position_dodge(width = 0.7), width = 0.2, colour = "grey40"
  ) +
  scale_fill_manual(values = GROUP_COLORS, name = "Responder") +
  labs(
    title = "Response magnitude",
    subtitle = "Median |logFC| per protein; whisker to the 90th percentile",
    x = NULL, y = "|logFC|"
  ) +
  FIG_THEME

# Panel B: HR-LR concordance per condition with Fisher-z 95% CI. Both weak; the CIs
# sit well below any moderate-agreement line, so the responders adapt largely apart.
pb <- ggplot(concordance, aes(phase, rho)) +
  geom_hline(yintercept = 0, colour = "grey80") +
  geom_pointrange(aes(ymin = rho_lo, ymax = rho_hi), colour = GROUP_COLORS[["HR"]]) +
  geom_text(
    aes(label = sprintf("\u03c1 = %.2f", rho)),
    vjust = -1.2, size = FIG_GEOM_TEXT, fontface = "bold"
  ) +
  scale_y_continuous(limits = c(-0.1, 0.6)) +
  labs(
    title = "HR-LR concordance",
    subtitle = "Spearman \u03c1 of per-protein logFC, Fisher-z 95% CI",
    x = NULL, y = "\u03c1 (HR vs LR)"
  ) +
  FIG_THEME

fig <- pa + pb +
  plot_annotation(
    title = "F00 pilot: proteome response is modest and the responders move largely apart",
    subtitle = "Training (T1\u2192T2) and acute (T2\u2192T3), HR vs LR, non-imputed limma results",
    theme = theme(
      plot.title = element_text(face = "bold", size = FIG_TITLE_SIZE),
      plot.subtitle = element_text(face = "italic", size = FIG_SUBTITLE_SIZE, colour = "grey30")
    )
  )

save_panel(fig, file.path(rpt, "magnitude_concordance"), 220, 110)

summary_tbl <- concordance |>
  left_join(
    magnitude |>
      pivot_wider(names_from = group, values_from = c(med_abs, p90_abs)),
    by = "phase"
  )
write_csv(summary_tbl, file.path(dat, "magnitude_concordance.csv"))

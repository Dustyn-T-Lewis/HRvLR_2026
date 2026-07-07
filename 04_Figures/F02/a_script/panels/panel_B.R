# F02 Panel B: Effect Size Distribution (logFC histograms)
# Signed logFC per contrast, coloured by the contrast code, density overlay,
# median |logFC| annotated. Faceted (stacked), YvO-style.

pacman::p_load(here, dplyr, tidyr, ggplot2)

if (!exists("meta")) source(here("04_Figures", "F02", "a_script", "HRvLR_F02_setup.R"))

PB_W <- 290
PB_H <- PANEL_H

display_contrasts <- setdiff(MAIN_CONTRASTS, c("Training_Interaction", "Acute_Interaction"))

lfc_long <- dep_df |>
  select(gene, starts_with("logFC_")) |>
  pivot_longer(starts_with("logFC_"), names_to = "contrast", values_to = "logFC") |>
  mutate(contrast = sub("logFC_", "", contrast)) |>
  filter(contrast %in% display_contrasts, !is.na(logFC)) |>
  mutate(contrast = factor(contrast, levels = display_contrasts))

lfc_stats <- lfc_long |>
  group_by(contrast) |>
  summarise(med_abs_lfc = median(abs(logFC)), .groups = "drop") |>
  mutate(annotation = sprintf("%.2f", med_abs_lfc))

x_lim <- as.numeric(quantile(abs(lfc_long$logFC), 0.99))
n_bins <- 50
binwidth <- 2 * x_lim / n_bins

pB <- ggplot(lfc_long, aes(logFC, fill = contrast)) +
  geom_histogram(binwidth = binwidth, color = "black", linewidth = 0.15, alpha = 0.9) +
  geom_density(aes(y = after_stat(count) * binwidth),
    alpha = 0.12, linewidth = 0.4, color = "grey20"
  ) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey50") +
  geom_text(
    data = lfc_stats, aes(x = -x_lim, y = Inf, label = annotation),
    inherit.aes = FALSE, hjust = 0, vjust = 1.4, size = FIG_GEOM_TEXT,
    color = "grey20", fontface = "bold"
  ) +
  facet_wrap(~contrast,
    nrow = 1, scales = "free_y",
    labeller = labeller(contrast = CTR_SHORT)
  ) +
  coord_cartesian(xlim = c(-x_lim, x_lim)) +
  scale_fill_manual(values = CONTRAST_COLORS) +
  labs(
    title = "Effect Size Distributions",
    subtitle = "logFC per contrast (5 contrasts) | value = median |logFC|",
    x = expression(log[2] ~ FC), y = "Count", tag = "B"
  ) +
  FIG_THEME +
  theme(legend.position = "none")

save_png(pB, file.path(RPT_DIR, "panels", "panel_b_effect_size"), PB_W, PB_H)
write.csv(lfc_stats, file.path(DAT_DIR, "audit_panel_B_lfc_stats.csv"), row.names = FALSE)
cat("F02 Panel B done.\n")

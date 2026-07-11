# F02 Panel D: effect-size histogram (logFC, down left / up right, counts on y)
# Pooled logFC histogram across the five main-effect contrasts with per-contrast
# density curves overlaid, so overall magnitude and per-contrast shape read in one
# compact panel. Title drawn on the composite.

pacman::p_load(here, dplyr, tidyr, ggplot2)

if (!exists("meta")) source(here("04_Figures", "F02", "a_script", "HRvLR_F02_setup.R"))

PD_W <- 120
PD_H <- 95

display_contrasts <- setdiff(MAIN_CONTRASTS, c("Training_Interaction", "Acute_Interaction"))

lfc_long <- dep_df |>
  select(gene, starts_with("logFC_")) |>
  pivot_longer(starts_with("logFC_"), names_to = "contrast", values_to = "logFC") |>
  mutate(contrast = sub("logFC_", "", contrast)) |>
  filter(contrast %in% display_contrasts, !is.na(logFC)) |>
  mutate(contrast = factor(contrast, levels = display_contrasts))

x_lim <- as.numeric(quantile(abs(lfc_long$logFC), 0.99))
binwidth <- 2 * x_lim / 45

pD <- ggplot(lfc_long, aes(logFC)) +
  geom_histogram(binwidth = binwidth, fill = "grey80", color = "white", linewidth = 0.1) +
  geom_density(aes(y = after_stat(count) * binwidth, color = contrast),
    linewidth = 0.5, key_glyph = "path"
  ) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey40") +
  scale_color_manual(values = CONTRAST_COLORS, labels = CTR_SHORT, name = NULL) +
  coord_cartesian(xlim = c(-x_lim, x_lim)) +
  labs(x = expression(log[2] ~ FC), y = "Proteins") +
  FIG_THEME +
  theme(
    legend.position = "bottom",
    legend.key.height = unit(2, "mm"),
    legend.text = element_text(size = FIG_LEGEND_TEXT - 1)
  ) +
  guides(color = guide_legend(nrow = 2, override.aes = list(linewidth = 1)))

save_png(pD, file.path(RPT_DIR, "panels", "panel_d_effectsize"), PD_W, PD_H)
write.csv(
  lfc_long |> group_by(contrast) |> summarise(med_abs_lfc = median(abs(logFC)), .groups = "drop"),
  file.path(DAT_DIR, "audit_panel_D_effectsize.csv"),
  row.names = FALSE
)
cat("F02 Panel D done.\n")

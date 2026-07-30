# F02 Panel E: effect-size distribution (per-contrast logFC density ridges)
# One horizontal ridge per contrast, contrasts on the y-axis in the same order as the
# DEP and pathway panels so the rows align across the bottom row. logFC left/right,
# median |logFC| annotated. Title drawn on the composite.

pacman::p_load(here, dplyr, tidyr, ggplot2, ggridges)

if (!exists("meta")) source(here("05_Figures", "F02_proteome", "a_script", "setup.R"))

PD_W <- 120
PD_H <- 110

display_contrasts <- setdiff(MAIN_CONTRASTS, c("Training_Interaction", "Acute_Interaction"))

lfc_long <- dep_df |>
  select(gene, starts_with("logFC_")) |>
  pivot_longer(starts_with("logFC_"), names_to = "contrast", values_to = "logFC") |>
  mutate(contrast = sub("logFC_", "", contrast)) |>
  filter(contrast %in% display_contrasts, !is.na(logFC)) |>
  mutate(contrast = factor(contrast, levels = rev(display_contrasts)))

x_lim <- as.numeric(quantile(abs(lfc_long$logFC), 0.99))
lfc_stats <- lfc_long |>
  group_by(contrast) |>
  summarise(med_abs_lfc = median(abs(logFC)), .groups = "drop")

pE <- ggplot(lfc_long, aes(logFC, contrast, fill = contrast)) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey50") +
  ggridges::geom_density_ridges(scale = 1.5, linewidth = 0.3, color = "grey25", alpha = 0.9) +
  geom_text(
    data = lfc_stats, aes(x = x_lim, y = contrast, label = sprintf("%.2f", med_abs_lfc)),
    inherit.aes = FALSE, hjust = 1, vjust = -0.4, size = FIG_GEOM_TEXT,
    fontface = "bold", color = "grey20"
  ) +
  scale_fill_manual(values = CONTRAST_COLORS, guide = "none") +
  scale_y_discrete(labels = CTR_SHORT, expand = expansion(add = c(0.1, 1.2))) +
  coord_cartesian(xlim = c(-x_lim, x_lim)) +
  labs(x = expression(log[2] ~ FC), y = NULL) +
  FIG_THEME +
  theme(axis.text.y = element_text(face = "bold"))

save_png(pE, file.path(RPT_DIR, "panels", "panel_e_effectsize"), PD_W, PD_H)
F02_AUDIT[["panel_E_effectsize"]] <- lfc_stats
cat("F02 Panel E done.\n")

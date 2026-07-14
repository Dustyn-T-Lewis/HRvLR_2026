# F02 DEPs-per-contrast panel, renders as composite tag D (diverging down/up, p and Pi counts)
# Down (left) / up (right) % of proteome; per direction two bars overlaid from zero:
# nominal p < 0.05 (light) and Pi < 0.05 (dark). Per-contrast background bands carry the
# contrast code; Pi counts sit at the bar tips.

pacman::p_load(here, dplyr, tidyr, tibble, ggplot2)

if (!exists("meta")) source(here("04_Figures", "F02_proteome", "a_script", "setup.R"))

PF_W <- 120
PF_H <- 110
n_total <- nrow(dep_df)

display_contrasts <- setdiff(MAIN_CONTRASTS, c("Training_Interaction", "Acute_Interaction"))

frac_df <- lapply(display_contrasts, function(ct) {
  p <- dep_df[[paste0("P.Value_", ct)]]
  lfc <- dep_df[[paste0("logFC_", ct)]]
  pi <- dep_df[[paste0("sig_pi_", ct)]]
  tibble(
    contrast = ct,
    key = c("Down p", "Down Pi", "Up p", "Up Pi"),
    direction = c("Down", "Down", "Up", "Up"),
    threshold = c("p", "Pi", "p", "Pi"),
    n = c(
      sum(p < 0.05 & lfc < 0, na.rm = TRUE), sum(pi == -1, na.rm = TRUE),
      sum(p < 0.05 & lfc > 0, na.rm = TRUE), sum(pi == 1, na.rm = TRUE)
    )
  )
}) |>
  bind_rows() |>
  mutate(
    pct = 100 * n / n_total,
    signed = if_else(direction == "Down", -pct, pct),
    contrast = factor(contrast, levels = rev(display_contrasts)),
    key = factor(key, levels = c("Down p", "Up p", "Down Pi", "Up Pi"))
  ) |>
  arrange(contrast, key)

DIR_FILL <- c(
  "Down p" = scales::alpha(DIR_COLORS[["Down"]], 0.40), "Down Pi" = DIR_COLORS[["Down"]],
  "Up p" = scales::alpha(DIR_COLORS[["Up"]], 0.40), "Up Pi" = DIR_COLORS[["Up"]]
)

band_levels <- levels(frac_df$contrast)
band_df <- tibble(
  xmin = seq_along(band_levels) - 0.5, xmax = seq_along(band_levels) + 0.5,
  band = CONTRAST_COLORS[band_levels]
)
pi_lab <- frac_df |> filter(threshold == "Pi", n > 0)
chance_pct <- 100 * 0.05 / 2

pE <- ggplot(frac_df, aes(contrast, signed, fill = key)) +
  geom_rect(
    data = band_df, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE, fill = scales::alpha(band_df$band, 0.14),
    color = "grey75", linewidth = 0.2
  ) +
  geom_col(position = "identity", width = 0.7, color = "white", linewidth = 0.3) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "grey30") +
  geom_hline(
    yintercept = c(-chance_pct, chance_pct), linetype = "dotted",
    linewidth = 0.3, color = "grey45"
  ) +
  geom_text(
    data = pi_lab, aes(contrast, signed, label = n),
    inherit.aes = FALSE, hjust = ifelse(pi_lab$direction == "Down", 1.25, -0.25),
    size = FIG_GEOM_TEXT, fontface = "bold", color = "grey15"
  ) +
  scale_fill_manual(
    values = DIR_FILL,
    breaks = c("Down Pi", "Up Pi"), labels = c("Down", "Up"), name = NULL
  ) +
  scale_x_discrete(labels = CTR_SHORT) +
  scale_y_continuous(labels = function(v) abs(v)) +
  coord_flip(clip = "off") +
  labs(x = NULL, y = "% of proteome") +
  FIG_THEME +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(face = "bold"),
    plot.margin = margin(6, 6, 4, 4)
  )

save_png(pE, file.path(RPT_DIR, "panels", "panel_e_dep_counts"), PF_W, PF_H)
F02_AUDIT[["panel_E_dep_counts"]] <- frac_df
cat("F02 Panel E done.\n")

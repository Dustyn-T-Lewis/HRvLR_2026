# F02 Panel B, continuous tree: DEPs per contrast (diverging down/up, p
# and Pi counts). Mirrors
# categorical/F02_proteome/a_script/panels/panel_d_dep_counts.R on the two
# pooled contrasts instead of nine, using POOLED_CONTRAST_COLORS/
# POOLED_CTR_SHORT (shared_style.R) in place of the group-coded palette.

pacman::p_load(here, dplyr, tidyr, tibble, ggplot2)

if (!exists("meta")) {
  source(here(
    "03_Analysis", "continuous", "F02_proteome", "a_script", "setup.R"
  ))
}

PF_W <- 120
PF_H <- 70
n_total <- nrow(dep_df)

frac_df <- lapply(MAIN_CONTRASTS, function(ct) {
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
    contrast = factor(contrast, levels = rev(MAIN_CONTRASTS)),
    key = factor(key, levels = c("Down p", "Up p", "Down Pi", "Up Pi"))
  ) |>
  arrange(contrast, key)

DIR_FILL <- c(
  "Down p" = scales::alpha(DIR_COLORS[["Down"]], 0.40),
  "Down Pi" = DIR_COLORS[["Down"]],
  "Up p" = scales::alpha(DIR_COLORS[["Up"]], 0.40), "Up Pi" = DIR_COLORS[["Up"]]
)

band_levels <- levels(frac_df$contrast)
band_df <- tibble(
  xmin = seq_along(band_levels) - 0.5, xmax = seq_along(band_levels) + 0.5,
  band = POOLED_CONTRAST_COLORS[band_levels]
)
pi_lab <- frac_df |> filter(threshold == "Pi", n > 0)
chance_pct <- 100 * 0.05 / 2

pB <- ggplot(frac_df, aes(contrast, signed, fill = key)) +
  geom_rect(
    data = band_df, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE, fill = scales::alpha(band_df$band, 0.14),
    color = "grey75", linewidth = 0.2
  ) +
  geom_col(
    position = "identity", width = 0.7, color = "white", linewidth = 0.3
  ) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "grey30") +
  geom_hline(
    yintercept = c(-chance_pct, chance_pct), linetype = "dotted",
    linewidth = 0.3, color = "grey45"
  ) +
  geom_text(
    data = pi_lab, aes(contrast, signed, label = n),
    inherit.aes = FALSE,
    hjust = ifelse(pi_lab$direction == "Down", 1.25, -0.25),
    size = FIG_GEOM_TEXT, fontface = "bold", color = "grey15"
  ) +
  scale_fill_manual(
    values = DIR_FILL,
    breaks = c("Down Pi", "Up Pi"), labels = c("Down", "Up"), name = NULL
  ) +
  scale_x_discrete(labels = POOLED_CTR_SHORT) +
  scale_y_continuous(labels = function(v) abs(v)) +
  coord_flip(clip = "off") +
  labs(x = NULL, y = "% of proteome") +
  FIG_THEME +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(face = "bold"),
    plot.margin = margin(6, 6, 4, 4)
  )

save_png(pB, file.path(RPT_DIR, "panels", "panel_b_dep_counts"), PF_W, PF_H)
F02_AUDIT[["panel_B_dep_counts"]] <- frac_df
cat("F02 (continuous) Panel B done.\n")

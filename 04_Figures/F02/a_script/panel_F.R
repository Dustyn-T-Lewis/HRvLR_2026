# F02 Panel F: DEPs per contrast (diverging up/down, nested p / Pi) ----
# Down (left) / up (right) % of proteome; each direction nested: nominal p < 0.05
# (light) then Pi < 0.05 (dark). Per-contrast background bands carry the contrast code.

setwd(rprojroot::find_rstudio_root_file())
if (!exists("meta")) source("04_Figures/F02/a_script/HRvLR_F02_setup.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})

PF_W <- 130
PF_H <- 120
n_total <- nrow(dep_df)

frac_df <- lapply(MAIN_CONTRASTS, function(c) {
  p <- dep_df[[paste0("P.Value_", c)]]
  lfc <- dep_df[[paste0("logFC_", c)]]
  pi <- dep_df[[paste0("sig_pi_", c)]]
  tibble(
    contrast = c,
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
    key = factor(key, levels = c("Down p", "Down Pi", "Up p", "Up Pi"))
  )

DIR_FILL <- c(
  "Down p" = scales::alpha(DIR_COLORS[["Down"]], 0.40), "Down Pi" = DIR_COLORS[["Down"]],
  "Up p" = scales::alpha(DIR_COLORS[["Up"]], 0.40), "Up Pi" = DIR_COLORS[["Up"]]
)

band_levels <- levels(frac_df$contrast)
band_df <- tibble(
  xmin = seq_along(band_levels) - 0.5, xmax = seq_along(band_levels) + 0.5,
  band = CONTRAST_COLORS[band_levels]
)
# Pi counts at the bar tips
pi_lab <- frac_df |> filter(threshold == "Pi", n > 0)

pF <- ggplot(frac_df, aes(contrast, signed, fill = key)) +
  geom_rect(
    data = band_df, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE, fill = scales::alpha(band_df$band, 0.14),
    color = "grey75", linewidth = 0.2
  ) +
  geom_col(position = "identity", width = 0.7, color = "black", linewidth = 0.2) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_text(
    data = pi_lab, aes(contrast, signed, label = n),
    inherit.aes = FALSE, hjust = ifelse(pi_lab$direction == "Down", 1.2, -0.2),
    size = 2.4, fontface = "bold"
  ) +
  scale_fill_manual(
    values = DIR_FILL,
    breaks = c("Down Pi", "Up Pi"), labels = c("Down", "Up"), name = NULL
  ) +
  scale_x_discrete(labels = CTR_SHORT) +
  coord_flip(clip = "off") +
  labs(
    title = "DEPs per Contrast",
    subtitle = sprintf(
      "%s proteins | down / up, nested p < 0.05 (light) + Pi < 0.05 (dark)",
      format(n_total, big.mark = ",")
    ),
    x = NULL, y = "% of proteome", tag = "F"
  ) +
  FIG_THEME +
  theme(legend.position = "bottom", axis.text.y = element_text(face = "bold"))

save_panel(pF, file.path(RPT_DIR, "panel_f_dep_counts"), PF_W, PF_H)
write.csv(frac_df, file.path(DAT_DIR, "audit_panel_F_dep_counts.csv"), row.names = FALSE)
cat("F02 Panel F done.\n")

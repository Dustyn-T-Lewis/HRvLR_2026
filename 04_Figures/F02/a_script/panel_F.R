# F02 Panel F: DEPs per contrast (diverging % of proteome, Pi < 0.05) ----
# Mirrors the YvO DEPs-per-contrast bar with per-contrast colour bands.

setwd(rprojroot::find_rstudio_root_file())
if (!exists("meta")) source("04_Figures/F02/a_script/HRvLR_F02_setup.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

PF_W <- 150
PF_H <- 120

n_total <- nrow(dep_df)

dep_counts <- lapply(MAIN_CONTRASTS, function(c) {
  sig <- dep_df[[paste0("sig_pi_", c)]]
  fdr <- dep_df[[paste0("adj.P.Val_", c)]]
  tibble(
    contrast = c,
    up = sum(sig == 1, na.rm = TRUE),
    down = sum(sig == -1, na.rm = TRUE),
    fdr = sum(fdr < 0.10, na.rm = TRUE)
  )
}) |>
  bind_rows() |>
  mutate(contrast = factor(contrast, levels = rev(MAIN_CONTRASTS)))

bar_df <- dep_counts |>
  pivot_longer(c(up, down), names_to = "direction", values_to = "n") |>
  mutate(
    pct = 100 * n / n_total,
    signed = if_else(direction == "down", -pct, pct),
    direction = factor(if_else(direction == "up", "Up", "Down"), c("Up", "Down"))
  )

pF <- ggplot(bar_df, aes(contrast, signed, fill = direction)) +
  geom_col(width = 0.74, color = "black", linewidth = 0.25) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_text(
    aes(
      label = ifelse(n > 0, n, ""),
      y = signed + sign(signed) * 0.15
    ),
    size = 2.6, hjust = ifelse(bar_df$signed >= 0, 0, 1)
  ) +
  scale_fill_manual(values = c(Up = unname(DIR_COLORS["Up"]), Down = unname(DIR_COLORS["Down"]))) +
  scale_x_discrete(labels = CTR_SHORT) +
  coord_flip(clip = "off") +
  labs(
    title = "DEPs per Contrast",
    subtitle = sprintf("%s proteins | Pi < 0.05 (down / up)", format(n_total, big.mark = ",")),
    x = NULL, y = "% of proteome", fill = NULL, tag = "F"
  ) +
  FIG_THEME +
  theme(legend.position = "bottom", axis.text.y = element_text(face = "bold"))

save_panel(pF, file.path(RPT_DIR, "panel_f_dep_counts"), PF_W, PF_H)
write.csv(dep_counts, file.path(DAT_DIR, "audit_panel_F_dep_counts.csv"), row.names = FALSE)
cat("F02 Panel F done.\n")

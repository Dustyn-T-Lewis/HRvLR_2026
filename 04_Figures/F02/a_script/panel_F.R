# F02 Panel F: DEPs per contrast (nested % of proteome) ----
# Per contrast, nested significance fractions: nominal p < 0.05 (lightest),
# FDR < 0.10, then Pi < 0.05 (darkest), shaded within the contrast colour.
# Mirrors the YvO DEPs-per-contrast bar.

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

frac_df <- lapply(MAIN_CONTRASTS, function(c) {
  p <- dep_df[[paste0("P.Value_", c)]]
  fdr <- dep_df[[paste0("adj.P.Val_", c)]]
  pi <- dep_df[[paste0("sig_pi_", c)]]
  tibble(
    contrast = c,
    threshold = c("p < 0.05", "FDR < 0.10", "Pi < 0.05"),
    n = c(
      sum(p < 0.05, na.rm = TRUE), sum(fdr < 0.10, na.rm = TRUE),
      sum(pi != 0, na.rm = TRUE)
    )
  )
}) |>
  bind_rows() |>
  mutate(
    pct = 100 * n / n_total,
    contrast = factor(contrast, levels = rev(MAIN_CONTRASTS)),
    threshold = factor(threshold, levels = c("p < 0.05", "FDR < 0.10", "Pi < 0.05")),
    fill_key = paste(contrast, threshold, sep = "___")
  )

# nested fill: contrast hue at increasing opacity per threshold
alpha_by <- c("p < 0.05" = 0.30, "FDR < 0.10" = 0.60, "Pi < 0.05" = 1)
FRAC_FILL <- vapply(levels(frac_df$fill_key %||% factor(frac_df$fill_key)), function(x) NA_character_, character(1))
FRAC_FILL <- setNames(
  mapply(
    function(c, t) scales::alpha(CONTRAST_COLORS[[as.character(c)]], alpha_by[[as.character(t)]]),
    frac_df$contrast, frac_df$threshold
  ),
  frac_df$fill_key
)

thr_lab <- c("p < 0.05" = "p", "FDR < 0.10" = "FDR", "Pi < 0.05" = "Π")
lab_df <- frac_df |>
  group_by(contrast) |>
  arrange(desc(pct), .by_group = TRUE) |>
  mutate(
    inner = lead(pct, default = 0), seg = pct - inner, lab_x = (pct + inner) / 2,
    label = thr_lab[as.character(threshold)]
  ) |>
  ungroup() |>
  filter(seg > 0.25)

pF <- ggplot(frac_df, aes(contrast, pct, fill = fill_key)) +
  geom_col(position = "identity", width = 0.74, color = "black", linewidth = 0.25) +
  geom_text(
    data = lab_df, aes(contrast, lab_x, label = label),
    inherit.aes = FALSE, size = 2.3, fontface = "bold", color = "grey15"
  ) +
  scale_fill_manual(values = FRAC_FILL, guide = "none") +
  scale_x_discrete(labels = CTR_SHORT) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.04))) +
  coord_flip(clip = "off") +
  labs(
    title = "DEPs per Contrast",
    subtitle = sprintf(
      "%s proteins | nested p < 0.05 / FDR < 0.10 / Pi < 0.05",
      format(n_total, big.mark = ",")
    ),
    x = NULL, y = "% of proteome", tag = "F"
  ) +
  FIG_THEME +
  theme(axis.text.y = element_text(face = "bold"))

save_panel(pF, file.path(RPT_DIR, "panel_f_dep_counts"), PF_W, PF_H)
write.csv(frac_df, file.path(DAT_DIR, "audit_panel_F_dep_counts.csv"), row.names = FALSE)
cat("F02 Panel F done.\n")

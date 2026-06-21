# F02 Panel B: Effect Size Distribution (|logFC| Histograms) ----
# Faceted by 5 MAIN_CONTRASTS, median |logFC| annotation

setwd(rprojroot::find_rstudio_root_file())
if (!exists("meta")) source("04_Figures/F02/a_script/HRvLR_F02_setup.R")

library(dplyr)
library(tidyr)
library(ggplot2)

PB_W <- 140; PB_H <- 160
dir.create(RPT_DIR, recursive = TRUE, showWarnings = FALSE)

lfc_long <- dep_df |>
  select(gene, starts_with("logFC_")) |>
  pivot_longer(starts_with("logFC_"), names_to = "contrast", values_to = "logFC") |>
  mutate(contrast = sub("logFC_", "", contrast)) |>
  filter(contrast %in% MAIN_CONTRASTS, !is.na(logFC)) |>
  mutate(abs_lfc = abs(logFC),
         contrast = factor(contrast, levels = MAIN_CONTRASTS))

lfc_stats <- lfc_long |>
  group_by(contrast) |>
  summarise(
    med_lfc     = median(logFC),
    med_abs_lfc = median(abs_lfc),
    n_above_05  = sum(abs_lfc > 0.5),
    .groups = "drop"
  )

pB <- ggplot(lfc_long, aes(x = abs_lfc)) +
  geom_histogram(bins = 50, fill = "grey60", color = "white", linewidth = 0.2) +
  geom_vline(data = lfc_stats, aes(xintercept = med_abs_lfc),
             linetype = "dashed", color = "red", linewidth = 0.5) +
  geom_text(data = lfc_stats,
            aes(x = med_abs_lfc, y = Inf,
                label = sprintf("median = %.3f", med_abs_lfc)),
            vjust = 1.5, hjust = -0.1, size = 2.5, color = "red") +
  facet_wrap(~ contrast, ncol = 1, scales = "free_y",
             labeller = labeller(contrast = CTR_SHORT)) +
  labs(title = "Effect Size Distributions",
       subtitle = "|logFC| per contrast (5 main)",
       x = "|logFC|", y = "Count", tag = "B") +
  FIG_THEME +
  theme(strip.text = element_text(size = 8, face = "bold"))

save_panel(pB, file.path(RPT_DIR, "panel_B_logfc_density_MAIN"), PB_W, PB_H)
write.csv(lfc_stats, file.path(DAT_DIR, "audit_panel_B_lfc_stats.csv"), row.names = FALSE)
cat("F02 Panel B done.\n")

# F03 Supplementary S1.7: Pi-score Distributions ----
# Histogram + ranked scatter per MAIN_CONTRAST

setwd(rprojroot::find_rstudio_root_file())
if (!exists("meta")) source("04_Figures/F03/a_script/HRvLR_F03_setup.R")

library(dplyr)
library(tidyr)
library(stringr)
library(patchwork)

dir.create(file.path(RPT_DIR, "supp"), recursive = TRUE, showWarnings = FALSE)

pi_long <- dep_df |>
  select(gene, starts_with("pi_score_")) |>
  pivot_longer(starts_with("pi_score_"), names_to = "contrast", values_to = "pi_score") |>
  mutate(contrast = str_remove(contrast, "pi_score_")) |>
  filter(contrast %in% MAIN_CONTRASTS, !is.na(pi_score))

p_hist <- ggplot(pi_long, aes(x = pi_score)) +
  geom_histogram(bins = 50, fill = "grey60", color = "white", linewidth = 0.2) +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.4) +
  facet_wrap(~ contrast, scales = "free_y", ncol = 2,
             labeller = labeller(contrast = CTR_SHORT)) +
  labs(x = expression(bold(Pi*"-score")), y = "Count") +
  FIG_THEME

pi_ranked <- pi_long |>
  group_by(contrast) |> arrange(pi_score) |>
  mutate(rank = row_number()) |> ungroup()

n_sig <- pi_long |>
  group_by(contrast) |>
  summarise(n = sum(pi_score < 0.05), .groups = "drop")

p_rank <- ggplot(pi_ranked, aes(x = rank, y = pi_score)) +
  geom_point(size = 0.2, alpha = 0.4, color = "grey40") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.4) +
  geom_text(data = n_sig, aes(label = sprintf("n = %d", n)),
            x = Inf, y = 0.10, hjust = 1.1, vjust = 0, size = 2.5, color = "red") +
  facet_wrap(~ contrast, scales = "free_x", ncol = 2,
             labeller = labeller(contrast = CTR_SHORT)) +
  labs(x = "Protein rank", y = expression(bold(Pi*"-score"))) +
  FIG_THEME

s17 <- p_hist / p_rank +
  plot_annotation(
    title = expression(bold("S1.7  ") * Pi * bold("-score distributions")),
    theme = theme(plot.title = element_text(size = 10))
  )

save_panel(s17, file.path(RPT_DIR, "supp", "supp_S1_7_SUPP"), 180, 200)
cat("F03 Supp S1.7 done.\n")

# Supplementary S1.7 — Pi-score Distribution Plots ----

if (!exists("meta")) source("04_Figures/F1/a_script/HRvLR_F1_setup.R")

S17_W <- 200
S17_H <- 160

# ── Reshape pi-scores to long format ----

pi_cols <- paste0("pi_score_", CONTRASTS)
pi_long <- dep_df |>
  select(gene, all_of(pi_cols)) |>
  pivot_longer(cols = all_of(pi_cols),
               names_to = "contrast",
               names_prefix = "pi_score_",
               values_to = "pi_score") |>
  filter(!is.na(pi_score)) |>
  mutate(contrast = factor(contrast, levels = CONTRASTS))

# Count DEPs (pi_score < 0.05) per contrast
n_sig <- pi_long |>
  group_by(contrast) |>
  summarise(n_dep = sum(pi_score < 0.05), .groups = "drop")

# ── Panel 1: Faceted histogram ----

p_hist <- ggplot(pi_long, aes(x = pi_score)) +
  geom_histogram(bins = 50, fill = "grey60", color = "white", linewidth = 0.2) +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.5) +
  facet_wrap(~ contrast, ncol = 4, labeller = labeller(contrast = CTR_FACET)) +
  labs(x = expression(Pi * "-score"),
       y = "Count") +
  FIG_THEME

# ── Panel 2: Faceted rank plot ----

pi_ranked <- pi_long |>
  group_by(contrast) |>
  arrange(desc(pi_score), .by_group = TRUE) |>
  mutate(rank = row_number()) |>
  ungroup()

# y position for n_sig label: slightly above the threshold line
label_df <- n_sig |>
  mutate(label = paste0("n = ", n_dep))

p_rank <- ggplot(pi_ranked, aes(x = rank, y = pi_score)) +
  geom_point(color = "grey40", size = 0.3, alpha = 0.6) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.5) +
  geom_text(data = label_df,
            aes(x = Inf, y = 0.08, label = label),
            inherit.aes = FALSE,
            hjust = 1.1, vjust = 0, size = 3, color = "red", fontface = "bold") +
  facet_wrap(~ contrast, ncol = 4, labeller = labeller(contrast = CTR_FACET),
             scales = "free_x") +
  labs(x = "Rank",
       y = expression(Pi * "-score")) +
  FIG_THEME

# ── Compose ----

s17 <- p_hist / p_rank +
  plot_annotation(
    title = expression(bold("S1.7  ") * Pi * bold("-score distributions"))
  )

# ── Export ----

save_panel(s17, file.path(RPT_DIR, "supplementary", "supp_S1_7"),
           width = S17_W, height = S17_H)

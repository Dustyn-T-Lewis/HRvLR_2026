# F02 Panel E: responder divergence, protein level (HR vs LR logFC, interaction-gated)
# Per protein, the training/acute logFC in HR against the same in LR. On-diagonal =
# concordant response; off-diagonal = divergent. Proteins whose divergence is
# significant carry it directly from the interaction contrast (sig_pi on
# Training_Interaction / Acute_Interaction), the correct test of HR-vs-LR difference
# that set overlap cannot provide (Gelman & Stern 2006). Divergent proteins labelled.

pacman::p_load(here, dplyr, tidyr, tibble, ggplot2, ggrepel)

if (!exists("meta")) source(here("04_Figures", "F02", "a_script", "HRvLR_F02_setup.R"))

PE_W <- 120
PE_H <- 110

phase_levels <- c("Training (T2 - T1)", "Acute (T3 - T2)")

diverg_df <- lapply(c(Training = "Training", Acute = "Acute"), function(ph) {
  tibble(
    gene = dep_df$gene,
    hr = dep_df[[paste0("logFC_", ph, "_HR")]],
    lr = dep_df[[paste0("logFC_", ph, "_LR")]],
    sig_int = dep_df[[paste0("sig_pi_", ph, "_Interaction")]]
  ) |>
    filter(!is.na(gene), !is.na(hr), !is.na(lr)) |>
    mutate(
      phase = ph,
      divergent = !is.na(sig_int) & sig_int != 0,
      gap = abs(hr - lr)
    )
}) |>
  bind_rows() |>
  mutate(phase = factor(recode(phase, Training = phase_levels[1], Acute = phase_levels[2]), levels = phase_levels))

label_df <- diverg_df |>
  filter(divergent) |>
  group_by(phase) |>
  slice_max(gap, n = 5, with_ties = FALSE) |>
  ungroup()

count_df <- diverg_df |>
  group_by(phase) |>
  summarise(n_div = sum(divergent), .groups = "drop") |>
  mutate(label = sprintf("%d divergent (Pi)", n_div))

lim <- as.numeric(quantile(abs(c(diverg_df$hr, diverg_df$lr)), 0.995))
INT_COLOR <- unname(CONTRAST_COLORS[["Acute_Interaction"]])

pC <- ggplot(diverg_df, aes(hr, lr)) +
  facet_wrap(~phase) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.4) +
  geom_hline(yintercept = 0, linewidth = 0.2, color = "grey80") +
  geom_vline(xintercept = 0, linewidth = 0.2, color = "grey80") +
  geom_point(data = ~ filter(.x, !divergent), color = "grey75", alpha = 0.3, size = 0.7) +
  geom_point(data = ~ filter(.x, divergent), color = INT_COLOR, alpha = 0.9, size = 1.4) +
  geom_text_repel(
    data = label_df, aes(label = gene), size = FIG_GEOM_TEXT, fontface = "bold",
    color = INT_COLOR, min.segment.length = 0, segment.size = 0.2,
    box.padding = 0.5, point.padding = 0.2, max.overlaps = Inf, seed = 42
  ) +
  geom_label(
    data = count_df, aes(x = -lim, y = lim, label = label), inherit.aes = FALSE,
    hjust = 0, vjust = 1, size = FIG_GEOM_TEXT, fontface = "bold", color = INT_COLOR,
    fill = scales::alpha("white", 0.85), label.size = 0, label.padding = unit(1.5, "pt")
  ) +
  coord_cartesian(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  labs(
    x = expression(bold(log[2] * FC[HR])), y = expression(bold(log[2] * FC[LR]))
  ) +
  FIG_THEME

save_png(pC, file.path(RPT_DIR, "panels", "panel_c_divergence"), PE_W, PE_H)
write.csv(
  diverg_df |> filter(divergent) |> select(phase, gene, hr, lr, gap),
  file.path(DAT_DIR, "audit_panel_C_divergence.csv"),
  row.names = FALSE
)
cat("F02 Panel C done.\n")

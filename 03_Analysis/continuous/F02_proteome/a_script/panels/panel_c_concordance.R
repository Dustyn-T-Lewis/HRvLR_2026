# F02 Panel C, continuous tree: does a protein that moves during training
# also move acutely? Replaces
# categorical/F02_proteome/a_script/panels/panel_c_divergence.R's HR-vs-LR
# logFC scatter, which has no group-free form -- there is no second arm to
# plot against. The natural continuous-tree question in the same on/off
# diagonal idiom is Training logFC against Acute logFC for the same
# protein: on-diagonal means the acute response echoes the training
# response, off-diagonal means the two phases pull apart. Candidates are
# proteins pi-gated on both contrasts, exploratory (pi carries no FDR
# control), same caveat categorical's panel states.
#
# This panel drops categorical's group-heterogeneity supplement
# (PERMDISP + per-arm CV scatter) -- both are inherently HR-vs-LR
# comparisons with no continuous analog.

pacman::p_load(here, dplyr, tidyr, tibble, ggplot2, ggrepel)

if (!exists("meta")) {
  source(here(
    "03_Analysis", "continuous", "F02_proteome", "a_script", "setup.R"
  ))
}

PE_W <- 120
PE_H <- 110

conc_df <- tibble(
  gene = dep_df$gene,
  training = dep_df$logFC_Training, acute = dep_df$logFC_Acute,
  pi_training = dep_df$sig_pi_Training, pi_acute = dep_df$sig_pi_Acute
) |>
  filter(!is.na(gene), !is.na(training), !is.na(acute)) |>
  mutate(
    candidate = !is.na(pi_training) & pi_training != 0 &
      !is.na(pi_acute) & pi_acute != 0,
    gap = abs(training - acute)
  )

label_df <- conc_df |>
  filter(candidate) |>
  mutate(direction = if_else(training - acute > 0, "up", "down")) |>
  group_by(direction) |>
  slice_max(gap, n = 5, with_ties = FALSE) |>
  ungroup()

n_candidate <- sum(conc_df$candidate)
r <- cor(conc_df$training, conc_df$acute)
count_label <- sprintf("%d pi-gated both phases | r = %.2f", n_candidate, r)

lim <- max(abs(c(conc_df$training, conc_df$acute))) * 1.1
CONC_COLOR <- unname(POOLED_CONTRAST_COLORS[["Acute"]])

pC <- ggplot(conc_df, aes(training, acute)) +
  coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  geom_abline(
    slope = 1, intercept = 0, linetype = "dashed", color = "grey55",
    linewidth = 0.4
  ) +
  geom_hline(yintercept = 0, linewidth = 0.2, color = "grey80") +
  geom_vline(xintercept = 0, linewidth = 0.2, color = "grey80") +
  geom_point(
    data = ~ filter(.x, !candidate), color = "grey75", alpha = 0.3, size = 0.7
  ) +
  geom_point(
    data = ~ filter(.x, candidate), color = CONC_COLOR, alpha = 0.9, size = 1.4
  ) +
  geom_text_repel(
    data = label_df, aes(label = gene), size = FIG_GEOM_TEXT, fontface = "bold",
    color = CONC_COLOR, min.segment.length = 0, segment.size = 0.2,
    box.padding = 0.6, point.padding = 0.25, max.overlaps = Inf, seed = 42
  ) +
  annotate(
    "label",
    x = -Inf, y = Inf, label = count_label, hjust = -0.06, vjust = 1.4,
    size = FIG_GEOM_TEXT, fontface = "bold", color = CONC_COLOR,
    fill = scales::alpha("white", 0.85), label.size = 0,
    label.padding = unit(1.5, "pt")
  ) +
  labs(
    x = expression(bold(log[2] * FC["Training"])),
    y = expression(bold(log[2] * FC["Acute"]))
  ) +
  FIG_THEME

save_png(pC, file.path(RPT_DIR, "panels", "panel_c_concordance"), PE_W, PE_H)
F02_AUDIT[["panel_C_concordance"]] <- conc_df |>
  filter(candidate) |>
  select(gene, training, acute, gap)
cat("F02 (continuous) Panel C done.\n")

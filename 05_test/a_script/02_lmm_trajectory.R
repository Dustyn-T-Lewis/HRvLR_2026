# Per-pathway mixed model of the GSVA trajectory: does the pathway follow a
# different HR-vs-LR path over T1-T3? score ~ group * timepoint + (1|subject).
if (!exists("expr")) source(here::here("05_test", "a_script", "setup.R"))
scores <- read_csv(file.path(DAT, "gsva_scores.csv"), show_col_types = FALSE) |>
  tibble::column_to_rownames("pathway") |>
  as.matrix()

long <- scores |>
  as.data.frame() |>
  tibble::rownames_to_column("pathway") |>
  pivot_longer(-pathway, names_to = "sample", values_to = "score") |>
  left_join(meta, by = "sample")

fit_one <- function(pw_name) {
  d <- long[long$pathway == pw_name, ]
  out <- tryCatch(
    {
      m <- lmerTest::lmer(score ~ group * timepoint + (1 | subject), data = d)
      a <- anova(m)
      eff <- d |>
        group_by(group, timepoint) |>
        summarise(mu = mean(score), .groups = "drop") |>
        pivot_wider(names_from = c(group, timepoint), values_from = mu)
      hr_change <- (eff$HR_T3 - eff$HR_T1)
      lr_change <- (eff$LR_T3 - eff$LR_T1)
      tibble(
        pathway = pw_name,
        p_interaction = a["group:timepoint", "Pr(>F)"],
        p_group = a["group", "Pr(>F)"],
        p_time = a["timepoint", "Pr(>F)"],
        hr_minus_lr_slope = hr_change - lr_change
      )
    },
    error = function(e) NULL
  )
  out
}

res <- bind_rows(lapply(rownames(scores), fit_one))
res$fdr_interaction <- p.adjust(res$p_interaction, "BH")
res$fdr_group <- p.adjust(res$p_group, "BH")
write_csv(res, file.path(DAT, "lmm_trajectory.csv"))

traj <- long |>
  group_by(pathway, group, timepoint) |>
  summarise(mean = mean(score), se = sd(score) / sqrt(dplyr::n()), .groups = "drop")
top6 <- res$pathway[order(res$p_interaction)][1:6]

fig_traj <- ggplot(
  traj |> filter(pathway %in% top6) |> mutate(pathway = clean_pathway_name(pathway, 34)),
  aes(timepoint, mean, color = group, group = group)
) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.12, linewidth = 0.4) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2) +
  scale_color_manual(values = GROUP_COLORS, name = NULL) +
  facet_wrap(~pathway, scales = "free_y", ncol = 3) +
  labs(
    x = NULL, y = "GSVA score",
    title = "Pathway trajectories with the strongest HR-vs-LR divergence",
    subtitle = "Mixed model score ~ group x timepoint + (1|subject); mean +/- SE"
  ) +
  FIG_THEME
save_panel(fig_traj, file.path(RPT, "fig1_trajectories"), 220, 150)

rank_df <- res |>
  slice_min(p_interaction, n = 20, with_ties = FALSE) |>
  mutate(
    term = reorder(clean_pathway_name(pathway, 40), -log10(p_interaction)),
    bias = if_else(hr_minus_lr_slope > 0, "HR rises more", "LR rises more"),
    fdr_hit = fdr_interaction < 0.05
  )
fig_inter <- ggplot(rank_df, aes(-log10(p_interaction), term, color = bias)) +
  geom_segment(aes(x = 0, xend = -log10(p_interaction), yend = term), linewidth = 0.5) +
  geom_point(aes(shape = fdr_hit), size = 2.4) +
  geom_vline(xintercept = -log10(0.05), linetype = "dotted", color = "grey50") +
  scale_color_manual(values = c("HR rises more" = "#2166AC", "LR rises more" = "#B2182B"), name = NULL) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 8), guide = "none") +
  labs(
    x = "-log10 p (group x timepoint)", y = NULL,
    title = "Group x timepoint interaction per pathway",
    subtitle = sprintf(
      "%d pathways, %d at nominal p < 0.05, %d at FDR < 0.05",
      nrow(res), sum(res$p_interaction < 0.05, na.rm = TRUE),
      sum(res$fdr_interaction < 0.05, na.rm = TRUE)
    )
  ) +
  FIG_THEME +
  theme(axis.text.y = element_text(size = 6.5))
save_panel(fig_inter, file.path(RPT, "fig2_interaction"), 170, 150)

cat(sprintf(
  "trajectory LMM: %d pathways, interaction FDR<0.05: %d\n",
  nrow(res), sum(res$fdr_interaction < 0.05, na.rm = TRUE)
))

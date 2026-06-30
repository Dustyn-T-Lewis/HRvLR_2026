# Shared builder for the F01 Pre/Post x Group + Delta panels (B/C/D/E/H/I): paired
# bars per Group_Time with HR/LR within-group brackets, a Delta-by-Group bar with a
# between-group test, and the mixed-ANOVA subtitle. Returns the panel and a Delta
# audit table. Jitter is seeded so the panel is reproducible.
pacman::p_load(here, dplyr, tidyr, readr, patchwork, ggsignif, rstatix)

prepost_delta_panel <- function(value_col, title, y_lab, delta_lab, tag,
                                y_comma = TRUE) {
  meta <- read_csv(here("00_input", "HRvLR_meta.csv"), show_col_types = FALSE) |>
    filter(Timepoint %in% c("T1", "T2"), !is.na(.data[[value_col]])) |>
    mutate(
      subject_key = sub("_T[123]$", "", Col_ID),
      Group = factor(Group, levels = c("HR", "LR")),
      Timepoint = factor(Timepoint, levels = c("T1", "T2")),
      Group_Time = factor(Group_Time, levels = c("HR_T1", "HR_T2", "LR_T1", "LR_T2")),
      value = .data[[value_col]]
    )

  pheno_wide <- meta |>
    select(subject_key, Group, Timepoint, value) |>
    pivot_wider(names_from = Timepoint, values_from = value) |>
    mutate(delta = T2 - T1)

  anova_tbl <- as.data.frame(rstatix::anova_test(
    data = meta, dv = value, wid = subject_key, between = Group, within = Timepoint
  ))
  hr_wide <- filter(pheno_wide, Group == "HR")
  lr_wide <- filter(pheno_wide, Group == "LR")
  paired_hr <- t.test(hr_wide$T2, hr_wide$T1, paired = TRUE)
  paired_lr <- t.test(lr_wide$T2, lr_wide$T1, paired = TRUE)
  delta_test <- t.test(delta ~ Group, data = pheno_wide)

  anova_sub <- sprintf(
    "Group %s   Time %s   Interaction %s",
    fmt_p(anova_tbl$p[anova_tbl$Effect == "Group"]),
    fmt_p(anova_tbl$p[anova_tbl$Effect == "Timepoint"]),
    fmt_p(anova_tbl$p[anova_tbl$Effect == "Group:Timepoint"])
  )

  jitter <- position_jitter(width = 0.12, height = 0, seed = 42)
  y_labels <- if (y_comma) scales::label_comma() else waiver()
  fill_4 <- GROUP_FILL[c("HR_T1", "HR_T2", "LR_T1", "LR_T2")]
  sig_y_left <- bracket_pos(meta$value)

  left <- ggplot(meta, aes(x = Group_Time, y = value, fill = Group_Time)) +
    geom_bar(stat = "summary", fun = mean, width = 0.65, color = "black", linewidth = 0.3) +
    geom_errorbar(stat = "summary", fun.data = mean_se, width = 0.2, linewidth = 0.4) +
    geom_point(position = jitter, size = 1.2, alpha = 0.5, shape = 21, color = "black", stroke = 0.3) +
    geom_signif(
      comparisons = list(c("HR_T1", "HR_T2")), annotations = fmt_p(paired_hr$p.value),
      y_position = sig_y_left, textsize = KEY_TEXT, tip_length = 0.01
    ) +
    geom_signif(
      comparisons = list(c("LR_T1", "LR_T2")), annotations = fmt_p(paired_lr$p.value),
      y_position = sig_y_left, textsize = KEY_TEXT, tip_length = 0.01
    ) +
    annotate("text", x = 1.5, y = -Inf, label = "HR", vjust = 4.2, fontface = "bold", size = 3.2, color = "grey25") +
    annotate("text", x = 3.5, y = -Inf, label = "LR", vjust = 4.2, fontface = "bold", size = 3.2, color = "grey25") +
    scale_fill_manual(values = fill_4) +
    scale_x_discrete(labels = c(HR_T1 = "Pre", HR_T2 = "Post", LR_T1 = "Pre", LR_T2 = "Post")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.22)), labels = y_labels) +
    coord_cartesian(clip = "off") +
    labs(title = title, subtitle = anova_sub, y = y_lab, x = NULL, tag = tag) +
    FIG_THEME +
    theme(
      plot.subtitle = element_text(size = 7, color = "grey40", face = "italic"),
      plot.margin = margin(5, 5, 20, 5), legend.position = "none"
    )

  delta_colors <- c(HR = unname(GROUP_FILL["HR_T2"]), LR = unname(GROUP_FILL["LR_T2"]))
  right <- ggplot(pheno_wide, aes(x = Group, y = delta, fill = Group)) +
    geom_bar(stat = "summary", fun = mean, width = 0.55, color = "black", linewidth = 0.3) +
    geom_errorbar(stat = "summary", fun.data = mean_se, width = 0.15, linewidth = 0.4) +
    geom_point(position = jitter, size = 1.2, alpha = 0.5, shape = 21, color = "black", stroke = 0.3) +
    geom_signif(
      comparisons = list(c("HR", "LR")), annotations = fmt_p(delta_test$p.value),
      textsize = KEY_TEXT, tip_length = 0.02, y_position = bracket_pos(pheno_wide$delta)
    ) +
    scale_fill_manual(values = delta_colors) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.25)), labels = y_labels) +
    labs(y = delta_lab, x = NULL) +
    FIG_THEME +
    theme(legend.position = "none")

  audit <- pheno_wide |>
    group_by(Group) |>
    summarise(
      n = sum(!is.na(delta)),
      delta_mean = mean(delta, na.rm = TRUE),
      delta_sd = sd(delta, na.rm = TRUE),
      delta_sem = delta_sd / sqrt(n),
      .groups = "drop"
    ) |>
    mutate(delta_t_test_p = delta_test$p.value)

  list(plot = (left | right) + plot_layout(widths = c(0.65, 0.35)), audit = audit)
}

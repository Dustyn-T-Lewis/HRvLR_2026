# Panel A: accumulated volume load, the matched-work control. A non-significant
# t-test with the HR-minus-LR difference and its 95% CI, read as "no detectable
# difference" rather than proof the groups trained identically (equivalence at
# n=8/group is inconclusive; see the qmd).
if (!exists("meta")) source(here::here("04_Figures", "F01_phenotype", "a_script", "setup.R"))
pacman::p_load(ggplot2, ggsignif, scales)

# Group-mean bar, SEM whiskers, seeded jitter, the responder palette, and an x
# axis labelled with each group's n.
group_bar_layers <- function(group) {
  list(
    geom_bar(stat = "summary", fun = mean, width = 0.5, color = "black", linewidth = 0.3),
    geom_errorbar(stat = "summary", fun.data = mean_se, width = 0.2, linewidth = 0.4),
    geom_point(
      position = position_jitter(width = 0.14, height = 0, seed = 42),
      size = 1.5, alpha = 0.5, shape = 16, color = "grey30"
    ),
    scale_fill_manual(values = GROUP_COLORS),
    scale_x_discrete(labels = c(
      HR = sprintf("HR (n=%d)", sum(group == "HR")),
      LR = sprintf("LR (n=%d)", sum(group == "LR"))
    ))
  )
}

f01_subtitle <- theme(
  plot.subtitle = element_text(size = 6.5, face = "italic", color = "grey40"),
  legend.position = "none"
)

build_volume <- function(meta, tag = "A") {
  vl <- meta |>
    filter(Timepoint == "T3", !is.na(ACCUM_VL)) |>
    transmute(subject, Group, value = ACCUM_VL)

  tt <- t.test(value ~ Group, data = vl)
  diff_hl <- unname(tt$estimate[1] - tt$estimate[2])
  ci <- tt$conf.int

  p <- ggplot(vl, aes(Group, value, fill = Group)) +
    group_bar_layers(vl$Group) +
    geom_signif(
      comparisons = list(c("HR", "LR")), annotations = fmt_p(tt$p.value),
      textsize = KEY_TEXT, tip_length = 0.02, y_position = bracket_pos(vl$value)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18)), labels = scales::label_comma()) +
    coord_cartesian(clip = "off") +
    labs(
      title = "Training work performed",
      subtitle = sprintf(
        "No detectable difference (p = %.2f; HR-LR 95%% CI [%s, %s] kg)",
        tt$p.value, scales::label_comma()(round(ci[1])), scales::label_comma()(round(ci[2]))
      ),
      y = "Volume load (kg)", x = NULL, tag = tag
    ) +
    FIG_THEME +
    f01_subtitle

  audit <- vl |>
    group_by(Group) |>
    summarise(n = n(), mean = mean(value), sd = sd(value), sem = sd / sqrt(n), .groups = "drop") |>
    mutate(diff_hr_lr = diff_hl, ci_lo = ci[1], ci_hi = ci[2], p = tt$p.value)
  list(plot = p, audit = audit)
}

a_volume <- build_volume(meta, tag = "A")
save_panel(a_volume$plot, file.path(F01_RPT, "panels", "panel_a_volume"), 150, 100)
F01_PANELS[["volume"]] <- a_volume$plot
F01_AUDIT[["volume_load"]] <- a_volume$audit

# Builders for the F01 phenotype figure. Each returns a ggplot plus the tidy
# table behind it. The divergence panels run on one linear mixed model per
# outcome; the volume control and the responder continuum need no model.
pacman::p_load(
  here, dplyr, tidyr, readr, tibble, purrr, forcats,
  patchwork, ggplot2, ggsignif, effectsize, lmerTest, emmeans
)
source(here("04_Figures", "functions", "style.R"))

# Repeated-measures outcomes (T1 -> T2). The MyoVision_fCSA_* columns are left
# out: at 200-560 they sit 10-30x below plausible fibre CSA and correlate
# negatively with manual tracing, the signature of a fibre-count metric
# mislabeled as area (MyoVision reports CSA in thousands of um^2; Wen et al.
# 2018, doi:10.1152/japplphysiol.00762.2017).
MEASURES <- tibble::tribble(
  ~col, ~label, ~domain,
  "mCSA_Pre", "mCSA", "muscle",
  "fCSA_Mixed_Pre", "fCSA Mixed", "fibre",
  "fCSA_Type_I_Pre", "fCSA Type I", "fibre",
  "fCSA_Type_II_Pre", "fCSA Type II", "fibre",
  "X1RM_Leg_Pre", "1RM Leg Press", "strength",
  "X1RM._Ext_Pre", "1RM Extension", "strength"
)

f01_meta <- function() {
  read_csv(here("00_input", "HRvLR_meta.csv"), show_col_types = FALSE) |>
    mutate(
      subject = sub("_T[123]$", "", Col_ID),
      Group = factor(Group, levels = c("HR", "LR")),
      Timepoint = factor(Timepoint, levels = c("T1", "T2", "T3"))
    )
}

f01_composite_scores <- function(meta) {
  meta |>
    filter(Timepoint == "T2", !is.na(COMP.HYPERTROPHY)) |>
    transmute(subject, Group, value = readr::parse_number(COMP.HYPERTROPHY))
}

prepost_long <- function(meta, col) {
  meta |>
    filter(Timepoint %in% c("T1", "T2"), !is.na(.data[[col]])) |>
    transmute(
      subject, Group,
      Timepoint = factor(Timepoint, levels = c("T1", "T2")),
      value = .data[[col]]
    )
}

# The raw mixed model for one outcome, reused by the report tables and the
# diagnostics. Fitted on the response scale by REML with Satterthwaite df.
fit_lmm <- function(meta, col) {
  lmerTest::lmer(
    value ~ Group * Timepoint + (1 | subject),
    data = prepost_long(meta, col),
    REML = TRUE
  )
}

# The divergence estimate: each subject's T2-T1 change, then the HR-minus-LR
# standardized difference (Hedges g, >0 = HR gains more) with its 95% CI. At two
# timepoints this equals the Group x Timepoint mixed-model interaction (Twisk et
# al. 2018, doi:10.1016/j.conctc.2018.03.008); the hlm supplement shows the mixed
# model agrees. Fibre, muscle, and strength share the SD axis.
change_advantage <- function(meta, col) {
  w <- prepost_long(meta, col) |>
    pivot_wider(names_from = Timepoint, values_from = value) |>
    filter(!is.na(T1), !is.na(T2)) |>
    mutate(delta = T2 - T1)
  g <- effectsize::hedges_g(delta ~ Group, data = w)
  tibble(
    estimate = g$Hedges_g, ci_lo = g$CI_low, ci_hi = g$CI_high,
    p = t.test(delta ~ Group, data = w)$p.value
  )
}

# Every outcome's divergence, Holm-adjusted across the six-outcome family.
change_advantage_table <- function(meta) {
  purrr::map_dfr(seq_len(nrow(MEASURES)), function(i) {
    change_advantage(meta, MEASURES$col[i]) |>
      mutate(measure = MEASURES$label[i], domain = MEASURES$domain[i], .before = 1)
  }) |>
    mutate(p_holm = p.adjust(p, method = "holm"))
}

# The six raw models, named by outcome, for the easystats report and diagnostics.
lmm_models <- function(meta) {
  setNames(lapply(MEASURES$col, function(col) fit_lmm(meta, col)), MEASURES$label)
}

# Leave-one-subject influence on the divergence coefficient. A balanced 2x2 design
# gives every observation identical leverage, so check_model's leverage view is
# uninformative; this refit-based shift is the meaningful influence check.
lmm_loso_influence <- function(meta) {
  bind_rows(lapply(seq_len(nrow(MEASURES)), function(i) {
    d <- prepost_long(meta, MEASURES$col[i])
    d$z <- as.numeric(scale(d$value))
    beta <- function(x) lme4::fixef(lmerTest::lmer(z ~ Group * Timepoint + (1 | subject), data = x))[["GroupLR:TimepointT2"]]
    full <- beta(d)
    shift <- vapply(unique(d$subject), function(s) full - beta(d[d$subject != s, ]), numeric(1))
    tibble(
      measure = MEASURES$label[i],
      max_abs_shift = max(abs(shift)),
      dfbetas_guide = 2 / sqrt(dplyr::n_distinct(d$subject))
    )
  }))
}

# Shared layers for the two HR/LR bar panels: group-mean bar, SEM whiskers,
# seeded jitter, the responder palette, and an x axis labelled with each group's n.
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

# Panel A: accumulated volume load, the matched-work control. A non-significant
# t-test with the HR-minus-LR difference and its 95% CI, read as "no detectable
# difference" rather than proof the groups trained identically (equivalence at
# n=8/group is inconclusive; see the qmd).
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

# Panel B: every subject's composite score, ranked. Carries both the split (the
# HR/LR boundary) and its magnitude as one continuum, so the defining score needs
# only one panel. Circular by construction (the score defines the groups).
build_continuum <- function(meta, tag = "B") {
  comp <- f01_composite_scores(meta) |>
    arrange(value) |>
    mutate(subject = fct_inorder(subject))
  boundary <- (min(comp$value[comp$Group == "HR"]) + max(comp$value[comp$Group == "LR"])) / 2
  swc <- 0.2 * sd(comp$value)

  p <- ggplot(comp, aes(value, subject, color = Group)) +
    annotate("rect",
      xmin = boundary - swc, xmax = boundary + swc, ymin = -Inf, ymax = Inf,
      fill = "grey85", alpha = 0.5
    ) +
    geom_vline(xintercept = boundary, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_segment(aes(x = 0, xend = value, yend = subject), linewidth = 0.5) +
    geom_point(size = 2.4) +
    scale_color_manual(values = GROUP_COLORS, name = NULL) +
    scale_x_continuous(labels = function(x) paste0(x, "%")) +
    labs(
      title = "The responder split, as a continuum", tag = tag,
      subtitle = "Per-subject composite %change; dashed = HR/LR split, grey = small-difference zone",
      x = "Composite hypertrophy (%)", y = NULL
    ) +
    FIG_THEME +
    theme(
      plot.subtitle = element_text(size = 6.5, face = "italic", color = "grey40"),
      axis.text.y = element_text(size = 6.5),
      legend.position = c(0.85, 0.2),
      panel.grid.major.y = element_blank()
    )
  audit <- comp |>
    transmute(subject, Group, composite = value, boundary, swc)
  list(plot = p, audit = audit)
}

# Panel C: the divergence forest. Standardized HR-minus-LR change (Hedges g) per
# outcome with 95% CI and the Holm-adjusted p; muscle and fibre clear zero,
# strength straddles it.
build_forest <- function(meta, tag = "C") {
  tbl <- change_advantage_table(meta)
  top_down <- c(
    "fCSA Type I", "fCSA Type II", "fCSA Mixed", "mCSA",
    "1RM Extension", "1RM Leg Press"
  )
  df <- tbl |>
    mutate(
      measure = factor(measure, levels = rev(top_down)),
      p_label = vapply(p_holm, fmt_p, character(1))
    )

  p <- ggplot(df, aes(estimate, measure, color = domain)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.4) +
    geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi), orientation = "y", width = 0.25, linewidth = 0.5) +
    geom_point(size = 2.6) +
    geom_text(
      aes(x = ci_hi, label = p_label),
      hjust = -0.15, size = 2.3, color = "grey35", show.legend = FALSE
    ) +
    scale_color_manual(
      values = DOMAIN_COLORS, name = "Domain",
      labels = c(fibre = "Fibre", muscle = "Muscle", strength = "Strength")
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.22))) +
    labs(
      title = "Phenotype change by domain (change score)", tag = tag,
      subtitle = "HR-LR change in each measure's SD units (Hedges g); Holm-adjusted p",
      x = "HR - LR change (SD of the measure)", y = NULL
    ) +
    FIG_THEME +
    theme(
      plot.subtitle = element_text(size = 6.5, face = "italic", color = "grey40"),
      legend.position = c(0.88, 0.17),
      legend.background = element_rect(fill = scales::alpha("white", 0.7), color = NA),
      legend.title = element_text(size = 8, face = "bold"),
      panel.grid.major.y = element_blank()
    )
  list(plot = p, audit = tbl)
}

# Supplement: change score by group for every outcome, evenly faceted. Delta is
# where the groups part (LR fibre CSA falls below zero); the mixed-model
# interaction p rides in each strip.
build_secondary <- function(meta) {
  adv <- change_advantage_table(meta)
  rows <- lapply(seq_len(nrow(MEASURES)), function(i) {
    lab <- MEASURES$label[i]
    p_holm <- adv$p_holm[adv$measure == lab]
    prepost_long(meta, MEASURES$col[i]) |>
      pivot_wider(names_from = Timepoint, values_from = value) |>
      transmute(subject, Group,
        delta = T2 - T1,
        measure = lab, order = i,
        facet = sprintf("%s  (p %s)", lab, fmt_p(p_holm))
      )
  })
  df <- bind_rows(rows) |> mutate(facet = fct_reorder(facet, order))

  p <- ggplot(df, aes(Group, delta, fill = Group)) +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
    geom_bar(stat = "summary", fun = mean, width = 0.55, color = "black", linewidth = 0.3) +
    geom_errorbar(stat = "summary", fun.data = mean_se, width = 0.15, linewidth = 0.4) +
    geom_point(
      position = position_jitter(width = 0.14, height = 0, seed = 42),
      size = 1.1, alpha = 0.5, shape = 21, color = "black", stroke = 0.3
    ) +
    facet_wrap(~facet, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = GROUP_COLORS, name = NULL) +
    labs(
      title = "Change score by group across all outcomes",
      subtitle = "Δ = T2-T1; Holm-adjusted change-score p per facet",
      y = "Δ (Post - Pre)", x = NULL
    ) +
    FIG_THEME +
    theme(legend.position = "none", strip.text = element_text(size = 8))

  audit <- df |>
    group_by(measure, Group) |>
    summarise(
      n = sum(!is.na(delta)), delta_mean = mean(delta, na.rm = TRUE),
      delta_sem = sd(delta, na.rm = TRUE) / sqrt(n), .groups = "drop"
    )
  list(plot = p, audit = audit)
}

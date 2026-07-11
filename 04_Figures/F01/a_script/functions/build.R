# Builders for the F01 phenotype figure. Each returns a ggplot plus the tidy
# table behind it. The divergence panels run on one linear mixed model per
# outcome; the volume control and the responder continuum need no model.
pacman::p_load(
  here, dplyr, tidyr, readr, tibble, purrr, forcats,
  patchwork, ggplot2, ggsignif, lmerTest, emmeans
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

# The divergence estimate: an emmeans interaction contrast on the z-scored
# outcome, so it reads as the HR-minus-LR difference in T1->T2 change in SD units
# (>0 = HR gains more) and fibre, muscle, and strength share one axis.
lmm_change_advantage <- function(meta, col) {
  fit <- lmerTest::lmer(
    scale(value)[, 1] ~ Group * Timepoint + (1 | subject),
    data = prepost_long(meta, col),
    REML = TRUE
  )
  emmeans::emmeans(fit, ~ Group * Timepoint) |>
    emmeans::contrast(method = list(adv = c(-1, 1, 0, 0) - c(0, 0, -1, 1))) |>
    summary(infer = TRUE) |>
    as.data.frame() |>
    transmute(estimate, ci_lo = lower.CL, ci_hi = upper.CL, p = p.value)
}

lmm_interaction_table <- function(meta) {
  bind_rows(lapply(seq_len(nrow(MEASURES)), function(i) {
    lmm_change_advantage(meta, MEASURES$col[i]) |>
      mutate(measure = MEASURES$label[i], domain = MEASURES$domain[i], .before = 1)
  }))
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

# Panel A: accumulated volume load, the matched-work control. An effect size with
# its interval and a TOST equivalence result, rather than reading a
# non-significant t-test as proof the groups trained the same.
build_volume <- function(meta, tag = "A") {
  pacman::p_load(TOSTER)
  vl <- meta |>
    filter(Timepoint == "T3", !is.na(ACCUM_VL)) |>
    transmute(subject, Group, value = ACCUM_VL)

  pooled_sd <- sqrt(mean(tapply(vl$value, vl$Group, var)))
  tost <- TOSTER::t_TOST(
    formula = value ~ Group,
    data = vl,
    eqb = 0.5 * pooled_sd # equivalence bound: half a pooled SD
  )
  g <- tost$effsize["Hedges's g(av)", c("estimate", "lower.ci", "upper.ci")]
  tost_p <- max(tost$TOST[c("TOST Lower", "TOST Upper"), "p.value"])
  nhst_p <- t.test(value ~ Group, data = vl)$p.value

  p <- ggplot(vl, aes(Group, value, fill = Group)) +
    group_bar_layers(vl$Group) +
    geom_signif(
      comparisons = list(c("HR", "LR")), annotations = fmt_p(nhst_p),
      textsize = KEY_TEXT, tip_length = 0.02, y_position = bracket_pos(vl$value)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18)), labels = scales::label_comma()) +
    coord_cartesian(clip = "off") +
    labs(
      title = "Training work performed",
      subtitle = sprintf("No group difference (g = %.2f, TOST p = %.2f)", g[["estimate"]], tost_p),
      y = "Volume load (kg)", x = NULL, tag = tag
    ) +
    FIG_THEME +
    f01_subtitle

  audit <- vl |>
    group_by(Group) |>
    summarise(n = n(), mean = mean(value), sd = sd(value), sem = sd / sqrt(n), .groups = "drop") |>
    mutate(hedges_g = g[["estimate"]], g_ci_lo = g[["lower.ci"]], g_ci_hi = g[["upper.ci"]], tost_p = tost_p)
  list(plot = p, audit = audit)
}

# Panel B: the composite score that defines the groups. No test (circular);
# shows the size of the split.
build_composite <- function(meta, tag = "B") {
  comp <- f01_composite_scores(meta)
  p <- ggplot(comp, aes(Group, value, fill = Group)) +
    group_bar_layers(comp$Group) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.12)), labels = function(x) paste0(x, "%")) +
    labs(
      title = "The responder split",
      subtitle = "Mean z-scored %change, mCSA + fibre CSA; defines groups (circular)",
      y = "Composite hypertrophy (%)", x = NULL, tag = tag
    ) +
    FIG_THEME +
    f01_subtitle
  audit <- comp |>
    group_by(Group) |>
    summarise(n = n(), mean = mean(value), sd = sd(value), sem = sd / sqrt(n), .groups = "drop")
  list(plot = p, audit = audit)
}

# Panel C: every subject's composite score, ranked. Shows the dichotomy as the
# tails of one continuum and that the cut falls in an empty gap wider than the
# smallest-worthwhile-change band.
build_continuum <- function(meta, tag = "C") {
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
      title = "Responders are a continuum", tag = tag,
      subtitle = "Subjects ranked; grey = small-difference zone, empty at the split",
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

# Panel D: the mixed-model divergence forest. Standardized HR-minus-LR change
# advantage per outcome, ordered by effect; muscle and fibre clear zero, strength
# straddles it.
build_forest <- function(meta, tag = "D") {
  tbl <- lmm_interaction_table(meta)
  domain_colors <- c(fibre = "#6BAED6", muscle = unname(GROUP_COLORS["HR"]), strength = "grey50")
  top_down <- c(
    "fCSA Type I", "fCSA Type II", "fCSA Mixed", "mCSA",
    "1RM Extension", "1RM Leg Press"
  )
  df <- tbl |> mutate(measure = factor(measure, levels = rev(top_down)))

  p <- ggplot(df, aes(estimate, measure, color = domain)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.4) +
    geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi), orientation = "y", width = 0.25, linewidth = 0.5) +
    geom_point(size = 2.6) +
    scale_color_manual(
      values = domain_colors, name = "Domain",
      labels = c(fibre = "Fibre", muscle = "Muscle", strength = "Strength")
    ) +
    labs(
      title = "Phenotype change by domain (mixed model)", tag = tag,
      subtitle = "HR-LR change scaled to each measure's SD, so domains compare on one axis",
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
  rows <- lapply(seq_len(nrow(MEASURES)), function(i) {
    p_int <- lmm_change_advantage(meta, MEASURES$col[i])$p
    prepost_long(meta, MEASURES$col[i]) |>
      pivot_wider(names_from = Timepoint, values_from = value) |>
      transmute(subject, Group,
        delta = T2 - T1,
        measure = MEASURES$label[i], order = i,
        facet = sprintf("%s  (int. %s)", MEASURES$label[i], fmt_p(p_int))
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
      subtitle = "Δ = T2-T1; mixed-model interaction p per facet",
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

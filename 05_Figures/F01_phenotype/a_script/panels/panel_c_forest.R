# Panel C: the divergence forest. Standardized HR-minus-LR change (Hedges g) per
# outcome with 95% CI and the Holm-adjusted p; muscle and fibre clear zero,
# strength straddles it.
#
# This panel owns two supplements, because both exist to support the same
# divergence claim:
#   1. the per-outcome change scores the forest summarises
#   2. the mixed-model diagnostics showing the LMM agrees with the change score
if (!exists("meta")) source(here::here("05_Figures", "F01_phenotype", "a_script", "setup.R"))
pacman::p_load(ggplot2, forcats, scales, performance, see, patchwork)

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
      legend.position = "inside",
      legend.position.inside = c(0.88, 0.17),
      legend.background = element_rect(fill = scales::alpha("white", 0.7), color = NA),
      legend.title = element_text(size = 8, face = "bold"),
      panel.grid.major.y = element_blank()
    )
  list(plot = p, audit = tbl)
}

# Supplement 1: change score by group for every outcome, evenly faceted. Delta is
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
        facet = sprintf("%s  (%s)", lab, fmt_p(p_holm))
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

c_forest <- build_forest(meta, tag = "C")
save_panel(c_forest$plot, file.path(F01_RPT, "panels", "panel_c_forest"), 150, 120)
F01_PANELS[["forest"]] <- c_forest$plot
F01_AUDIT[["change_advantage"]] <- c_forest$audit

secondary <- build_secondary(meta)
rpt_supp <- file.path(F01_RPT, "supp")
dir.create(rpt_supp, recursive = TRUE, showWarnings = FALSE)
save_panel(secondary$plot, file.path(rpt_supp, "panel_c_supp_change_scores"), 230, 150)
F01_AUDIT[["secondary_delta"]] <- secondary$audit

# Supplement 2: one check_model exemplar (mCSA) and the per-model fit table. The
# six models share structure, so one figure shows the assumption pattern;
# compare_performance validates all six quantitatively.
models <- lmm_models(meta)

# check_model simulates from the fit for its posterior-predictive panel, so the
# figure is a different draw on every run unless the RNG is pinned.
set.seed(20260713)

ggsave(
  file.path(rpt_supp, "panel_c_supp_check_model.png"),
  plot(check_model(models[["mCSA"]], verbose = FALSE)),
  width = 280, height = 360, units = "mm", dpi = 150, bg = "white" # tall so y-titles clear the subtitles
)

F01_AUDIT[["lmm_fit_summary"]] <- as.data.frame(
  performance::compare_performance(models, metrics = c("AIC", "BIC", "R2", "ICC", "RMSE"))
) |>
  dplyr::transmute(
    measure = Name, AIC, BIC,
    r2_conditional = R2_conditional, r2_marginal = R2_marginal, icc = ICC, RMSE
  )

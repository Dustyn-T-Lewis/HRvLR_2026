# Panel builders for the prediction suite. Compact ROC / AUC and
# predicted-vs-observed / Q^2 figures at roughly F07 scale, styled through the
# shared 04_Figures style.R.

pacman::p_load(ggplot2, dplyr, tidyr, patchwork, scales)

SPACE_COLORS <- c(
  singscore = "#238B45", proteins = unname(GROUP_COLORS["HR"]),
  eigengenes = "#9E9AC8"
)
SPACE_LABELS <- c(
  singscore = "singscore", proteins = "proteins", eigengenes = "module ME"
)
MODEL_LABELS <- c(glmnet = "elastic net", spls = "sPLS")
MODEL_COLORS <- c(
  glmnet = unname(GROUP_COLORS["LR"]), spls = unname(GROUP_COLORS["HR"])
)

# Class arm: ROC small-multiples over feature space (top) above AUC bars with
# DeLong intervals, the chance reference at 0.5, permutation p and BH stars.
build_class_panel <- function(summ, roc_df, phase) {
  summ <- summ |>
    mutate(
      space = factor(space, levels = names(SPACE_COLORS)),
      lab = sprintf("%.2f\np=%.3f%s", estimate, perm_p, sig_stars(q_bh))
    )
  roc_df <- roc_df |> mutate(space = factor(space, levels = names(SPACE_COLORS)))

  p_roc <- ggplot(roc_df, aes(fpr, tpr, colour = model)) +
    geom_abline(
      slope = 1, intercept = 0, linetype = "dashed",
      colour = "grey70", linewidth = 0.3
    ) +
    geom_line(linewidth = 0.9) +
    facet_wrap(~space, nrow = 1, labeller = as_labeller(SPACE_LABELS)) +
    scale_colour_manual(
      values = MODEL_COLORS,
      labels = MODEL_LABELS, name = NULL
    ) +
    scale_x_continuous(breaks = c(0, 1)) +
    scale_y_continuous(breaks = c(0, 1)) +
    coord_equal(clip = "off") +
    labs(x = "1 - specificity", y = "sensitivity") +
    FIG_THEME +
    theme(legend.position = "top", panel.spacing = unit(6, "pt"))

  p_bar <- ggplot(summ, aes(space, estimate, fill = model)) +
    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
    geom_col(
      position = position_dodge(0.7), width = 0.62,
      colour = "black", linewidth = 0.25
    ) +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
      position = position_dodge(0.7), width = 0.2, linewidth = 0.35
    ) +
    geom_text(aes(y = pmax(ci_hi, estimate) + 0.05, label = lab),
      position = position_dodge(0.7), size = 2.3, lineheight = 0.85,
      fontface = "bold", colour = "grey15"
    ) +
    scale_fill_manual(
      values = MODEL_COLORS,
      labels = MODEL_LABELS, name = NULL
    ) +
    scale_x_discrete(labels = SPACE_LABELS) +
    scale_y_continuous(limits = c(0, 1.15), breaks = seq(0, 1, 0.25)) +
    labs(x = NULL, y = "LOSO AUC") +
    FIG_THEME +
    theme(legend.position = "none")

  (p_roc / p_bar) +
    plot_layout(heights = c(1, 1.1)) +
    plot_annotation(
      title = sprintf("Responder-class prediction — %s features", phase),
      subtitle = sprintf("Nested leave-one-subject-out AUC vs a %d-permutation null. Dashed line = chance; * BH q<.05.", N_PERM),
      theme = theme(
        plot.title = element_text(face = "bold", size = FIG_TITLE_SIZE),
        plot.subtitle = element_text(
          face = "italic", size = FIG_SUBTITLE_SIZE,
          colour = "grey30"
        )
      )
    )
}

# Continuous arm: Q^2 bars over feature space faceted by outcome (top) above the
# singscore predicted-vs-observed scatters (bottom). Null reference at Q^2 = 0.
build_cont_panel <- function(summ, preds, phase) {
  summ <- summ |>
    mutate(
      space = factor(space, levels = names(SPACE_COLORS)),
      q2_disp = pmin(pmax(q2, -0.6), 0.85),
      lab = sprintf("%.2f%s", q2, sig_stars(q_bh))
    )

  p_bar <- ggplot(summ, aes(space, q2_disp, fill = model)) +
    geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.4) +
    geom_col(
      position = position_dodge(0.7), width = 0.62,
      colour = "black", linewidth = 0.2
    ) +
    geom_text(aes(y = pmax(q2_disp, 0) + 0.04, label = lab),
      position = position_dodge(0.7), size = 2.0, fontface = "bold",
      colour = "grey15"
    ) +
    facet_wrap(~outcome, nrow = 1) +
    scale_fill_manual(
      values = MODEL_COLORS,
      labels = MODEL_LABELS, name = NULL
    ) +
    scale_x_discrete(
      labels = SPACE_LABELS,
      guide = guide_axis(n.dodge = 2)
    ) +
    scale_y_continuous(limits = c(-0.62, 0.92)) +
    labs(x = NULL, y = expression(LOSO ~ Q^2)) +
    FIG_THEME +
    theme(
      legend.position = "top", panel.spacing = unit(5, "pt"),
      axis.text.x = element_text(size = 6)
    )

  scatter_df <- preds |>
    filter(space == "singscore", model == "glmnet")
  p_sc <- ggplot(scatter_df, aes(y, pred)) +
    geom_abline(
      slope = 1, intercept = 0, linetype = "dashed",
      colour = "grey70", linewidth = 0.3
    ) +
    geom_point(colour = SPACE_COLORS[["singscore"]], size = 1.4, alpha = 0.8) +
    facet_wrap(~outcome, nrow = 1, scales = "free") +
    labs(x = "observed", y = "LOSO predicted") +
    FIG_THEME +
    theme(panel.spacing = unit(5, "pt"))

  (p_bar / p_sc) +
    plot_layout(heights = c(1.1, 1)) +
    plot_annotation(
      title = sprintf("Continuous-phenotype prediction — %s features", phase),
      subtitle = sprintf("Nested LOSO Q^2 vs a %d-permutation null; scatters are singscore + elastic net. Null reference at Q^2 = 0.", N_PERM),
      theme = theme(
        plot.title = element_text(face = "bold", size = FIG_TITLE_SIZE),
        plot.subtitle = element_text(
          face = "italic", size = FIG_SUBTITLE_SIZE,
          colour = "grey30"
        )
      )
    )
}

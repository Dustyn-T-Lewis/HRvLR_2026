# Panel builders for the F05 prediction leaves. Each leaf figure shows its
# observed metric against the permutation-null band for each feature space,
# with the honest-null framing on its face: a bar near chance with a
# non-significant permutation p is the expected, reported result.

pacman::p_load(ggplot2, dplyr, tidyr, patchwork, scales)

SPACE_LEVELS <- c("proteins", "singscore", "eigengenes")
SPACE_LABELS <- c(
  proteins = "proteins", singscore = "pathways (singscore)",
  eigengenes = "modules (ME)*"
)
SPACE_COLORS <- c(
  proteins = "#2166AC", singscore = "#238B45", eigengenes = "#9E9AC8"
)
MODEL_LABELS <- c(glmnet = "elastic net", spls = "sPLS-DA", pam = "PAM")

space_factor <- function(x) factor(x, levels = SPACE_LEVELS)

# Classification leaf: AUC bar per feature space with the DeLong interval, the
# permutation-null mean as a reference, and the top selected features by fold
# frequency beside it.
build_class_leaf_panel <- function(summ, roc_df, selection, contrast, model) {
  summ <- summ |>
    mutate(
      space = space_factor(space),
      lab = sprintf("AUC %.2f\np=%.3f", estimate, perm_p)
    )

  p_auc <- ggplot(summ, aes(space, estimate, fill = space)) +
    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
    geom_col(width = 0.62, colour = "black", linewidth = 0.25) +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
      width = 0.18, linewidth = 0.35
    ) +
    geom_point(aes(y = null_mean),
      shape = 18, size = 2.4, colour = "grey25"
    ) +
    geom_text(aes(y = pmax(ci_hi, estimate) + 0.06, label = lab),
      size = 2.4, lineheight = 0.85, fontface = "bold", colour = "grey15"
    ) +
    scale_fill_manual(values = SPACE_COLORS, guide = "none") +
    scale_x_discrete(labels = SPACE_LABELS) +
    scale_y_continuous(limits = c(0, 1.2), breaks = seq(0, 1, 0.25)) +
    labs(
      x = NULL, y = "LOSO AUC",
      caption = "◆ = permutation-null mean; dashed = chance"
    ) +
    FIG_THEME +
    theme(plot.caption = element_text(size = 7, colour = "grey40"))

  top_sel <- selection |>
    filter(space == "proteins") |>
    slice_max(freq, n = 8, with_ties = FALSE)
  p_sel <- if (nrow(top_sel)) {
    ggplot(top_sel, aes(freq, reorder(feature, freq))) +
      geom_col(fill = SPACE_COLORS[["proteins"]], width = 0.7) +
      scale_x_continuous(limits = c(0, 1)) +
      labs(x = "fold selection freq.", y = NULL, subtitle = "top proteins") +
      FIG_THEME
  } else {
    patchwork::plot_spacer()
  }

  (p_auc | p_sel) +
    plot_layout(widths = c(1.3, 1)) +
    plot_annotation(
      title = sprintf(
        "%s (%s) — %s", contrast, CONTRAST_ROLE[[contrast]],
        MODEL_LABELS[[model]]
      ),
      subtitle = sprintf(
        paste(
          "Nested LOSO vs %d-permutation null.",
          "*module eigengenes are cohort-relative (optimistic)."
        ),
        N_PERM
      ),
      theme = theme(
        plot.title = element_text(face = "bold", size = FIG_TITLE_SIZE),
        plot.subtitle = element_text(
          face = "italic", size = FIG_SUBTITLE_SIZE, colour = "grey30"
        )
      )
    )
}

# Continuous leaf: Q^2 per space faceted by outcome, null reference at 0.
build_cont_leaf_panel <- function(summ, preds, model) {
  summ <- summ |>
    mutate(
      space = space_factor(space),
      q2_disp = pmin(pmax(q2, -0.6), 0.85),
      lab = sprintf("%.2f\np=%.2f", q2, perm_p_q2)
    )

  ggplot(summ, aes(space, q2_disp, fill = space)) +
    geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.4) +
    geom_col(width = 0.62, colour = "black", linewidth = 0.2) +
    geom_text(aes(y = pmax(q2_disp, 0) + 0.05, label = lab),
      size = 2.0, lineheight = 0.85, fontface = "bold", colour = "grey15"
    ) +
    facet_wrap(~outcome, nrow = 1) +
    scale_fill_manual(values = SPACE_COLORS, guide = "none") +
    scale_x_discrete(labels = SPACE_LABELS, guide = guide_axis(n.dodge = 2)) +
    scale_y_continuous(limits = c(-0.62, 0.92)) +
    labs(
      x = NULL, y = expression(LOSO ~ Q^2),
      title = sprintf(
        "Baseline continuous prediction — %s", MODEL_LABELS[[model]]
      ),
      subtitle = sprintf(
        paste(
          "Q^2 <= 0 means no better than the outcome mean.",
          "Nested LOSO vs %d-permutation null."
        ),
        N_PERM
      )
    ) +
    FIG_THEME +
    theme(
      axis.text.x = element_text(size = 6),
      plot.title = element_text(face = "bold", size = FIG_TITLE_SIZE),
      plot.subtitle = element_text(
        face = "italic", size = FIG_SUBTITLE_SIZE, colour = "grey30"
      )
    )
}

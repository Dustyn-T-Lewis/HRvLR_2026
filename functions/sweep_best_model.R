# Best-performing cell per feature level, drawn against the screen it was
# selected from. Two rules can disagree -- the largest metric and the smallest
# permutation p -- so both winners are drawn and labelled with the rule that
# chose them. BH is monotone in p and adds no ranking, so q is reported as a
# pass/fail annotation rather than a third selection rule.

suppressWarnings(suppressMessages({
  library(here)
  source(here("functions", "sweep_composites.R"))
  source(here("functions", "sweep_cell_panel.R"))
  pacman::p_load(openxlsx, dplyr, ggplot2, patchwork, pROC)
}))

BEST_SPEC <- list(
  F05_classification = list(
    metric = "estimate", p = "perm_p", null_col = "auc", chance = 0.5,
    label = "AUC", xlab = "permutation-null AUC", panel = "roc"
  ),
  F06_prediction = list(
    metric = "q2", p = "perm_p_q2", null_col = "q2", chance = 0,
    label = "Q2", xlab = "permutation-null Q2", panel = "obs_pred"
  )
)

# One row per winning cell, tagged with every rule that selected it. A single
# cell winning on both rules yields one row, not two.
best_cells <- function(root, level) {
  spec <- BEST_SPEC[[root]]
  d <- root_cells(root) |>
    filter(.data$level == .env$level) |>
    mutate(q = p.adjust(.data[[spec$p]], "BH"))
  chosen <- c(
    metric = which.max(d[[spec$metric]]), p = which.min(d[[spec$p]])
  )
  lapply(unique(chosen), function(i) {
    d[i, ] |> mutate(
      rule = paste(names(chosen)[chosen == i], collapse = " + "),
      n_screened = nrow(d)
    )
  }) |>
    bind_rows()
}

best_null_panel <- function(root, cell) {
  spec <- BEST_SPEC[[root]]
  nulls <- leaf_sheet(
    root, cell$level, cell$config, cell_phenotype(cell), cell$model, "null"
  )
  best_null_panel_from(nulls[[spec$null_col]], cell[[spec$metric]], spec)
}

best_null_panel_from <- function(values, observed, spec) {
  null_panel(
    values, observed, spec$chance, spec$xlab,
    "observed (red) vs null (grey); dotted = the metric's own baseline"
  )
}

# The footer is the honest part of the figure: it names the screen the winner
# was drawn from and whether it survives correction across that screen.
best_footer <- function(root, cell) {
  spec <- BEST_SPEC[[root]]
  sprintf(
    paste(
      "%s = %.3f  ·  nominal p = %.4g  ·  BH q = %.3g (%s)",
      " ·  best of %d cells by %s"
    ),
    spec$label, cell[[spec$metric]], cell[[spec$p]], cell$q,
    if (cell$q < 0.05) "passes" else "fails", cell$n_screened, cell$rule
  )
}

best_cell_panel <- function(root, cell) {
  spec <- BEST_SPEC[[root]]
  main <- if (identical(spec$panel, "roc")) {
    roc_panel(root, cell)
  } else {
    obs_pred_panel(root, cell)
  }
  (main | best_null_panel(root, cell)) +
    plot_layout(widths = c(1, 1.2)) +
    plot_annotation(
      title = sprintf(
        "%s · %s · %s · %s", cell$level, cell$config,
        cell_phenotype(cell), cell$model
      ),
      subtitle = best_footer(root, cell),
      theme = leaf_title_theme()
    )
}

build_best_model_level <- function(root, level) {
  cells <- best_cells(root, level)
  out_dir <- file.path(sweep_root_dir(root), level, "best_model", "b_reports")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  for (i in seq_len(nrow(cells))) {
    save_panel(
      best_cell_panel(root, cells[i, ]),
      file.path(out_dir, paste0("best_", gsub(" \\+ ", "_", cells$rule[i]))),
      300, 110
    )
  }

  panels <- lapply(seq_len(nrow(cells)), function(i) {
    spec <- BEST_SPEC[[root]]
    cell <- cells[i, ]
    main <- if (identical(spec$panel, "roc")) {
      roc_panel(root, cell)
    } else {
      obs_pred_panel(root, cell)
    }
    (main | best_null_panel(root, cell)) +
      plot_annotation(subtitle = best_footer(root, cell))
  })

  composite <- wrap_plots(lapply(panels, patchwork::wrap_elements), ncol = 1) +
    plot_annotation(
      title = sprintf("%s -- best %s cell", root, level),
      subtitle = paste(
        "Winners by largest metric and by smallest permutation p. Where the",
        "two rules disagree both are drawn.\nEvery panel names the screen it",
        "was selected from; BH is taken across that screen."
      ),
      tag_levels = "A",
      theme = theme(
        plot.title = element_text(face = "bold", size = FIG_TITLE_SIZE + 2),
        plot.subtitle = element_text(
          face = "italic", size = FIG_SUBTITLE_SIZE, colour = "grey30"
        ),
        plot.tag = element_text(face = "bold", size = FIG_TAG_SIZE)
      )
    )
  save_panel(
    composite, file.path(out_dir, "best_model"), 300, 110 * nrow(cells) + 30
  )
  message(root, " ", level, ": ", nrow(cells), " winner(s) -> ", out_dir)
  cells
}

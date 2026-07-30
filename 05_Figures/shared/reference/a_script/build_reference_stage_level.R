# One reference set per stage x feature level, written beside that level in the
# reference tree. Each set holds two views in the YvO idiom:
#
#   heatmap  - the F06 view: features x outcomes, timepoint configs as column
#              families, filled by the stage's signed statistic.
#   detail   - the F07 view: the strongest cells shown as raw data, so a claim
#              can be checked against the observations behind it.

suppressWarnings(suppressMessages({
  library(here)
  source(here("functions", "sweep_composites.R"))
  source(here("functions", "sweep_drivers.R"))
  source(here("functions", "pred_features.R"))
  pacman::p_load(openxlsx, dplyr, tidyr, ggplot2, patchwork, stringr)
}))

REF_ROOT <- here("05_Figures", "shared", "reference")
SPARSE_MODELS <- c("lasso", "enet", "spls", "pam")
OUTCOMES <- c(ADAPT_OUTCOMES, "group_diff")
OUTCOME_SHORT <- c(
  comp_hypertrophy = "comp.hyp", d_fcsa_I = "fCSA I", d_fcsa_II = "fCSA II",
  d_mcsa = "mCSA", d_1rm_legpress = "1RM leg", d_1rm_ext = "1RM ext",
  group_diff = "HR-LR"
)
N_SHOW <- c(modules = 13, pathways = 22, proteins = 40)

strip_tp <- function(x) sub("@T[123]$", "", x)

cell_files <- function(stage, level, config = "*", model = "*") {
  Sys.glob(file.path(
    sweep_root_dir(stage), level, config, "*", model, "c_data", "results.xlsx"
  ))
}

row_labels <- function(features, level) {
  keys <- driver_keys(features, level)
  out <- if (level == "proteins") {
    features
  } else {
    paste0(keys$key, ": ", keys$description)
  }
  stats::setNames(make.unique(out), features)
}

# Prediction and classification: how often each feature is selected, pooled over
# the sparse models. Frequency is unsigned, so the scale runs 0 to 1.
select_matrix <- function(stage, level) {
  bind_rows(lapply(cell_files(stage, level), function(f) {
    sel <- read.xlsx(f, "selection")
    if (!nrow(sel)) {
      return(NULL)
    }
    if (!"outcome" %in% names(sel)) sel$outcome <- "group_diff"
    sel |>
      filter(.data$model %in% SPARSE_MODELS) |>
      mutate(
        feature = strip_tp(.data$feature),
        config = basename(dirname(dirname(dirname(dirname(f)))))
      )
  })) |>
    group_by(.data$feature, .data$outcome, .data$config) |>
    summarise(value = mean(.data$freq), .groups = "drop") |>
    mutate(p = NA_real_)
}

stage_matrix <- function(stage, level) {
  select_matrix(stage, level)
}

stage_heatmap <- function(stage, level) {
  d <- stage_matrix(stage, level)
  if (!nrow(d)) {
    return(NULL)
  }
  rank_by <- d |>
    group_by(.data$feature) |>
    summarise(s = -max(.data$value, na.rm = TRUE), .groups = "drop") |>
    arrange(.data$s)
  keep <- head(rank_by$feature, N_SHOW[[level]])
  labs_map <- row_labels(keep, level)

  d <- d |>
    filter(.data$feature %in% keep) |>
    mutate(
      row = factor(labs_map[.data$feature], levels = rev(labs_map)),
      config = factor(.data$config, levels = SWEEP_CONFIGS),
      outcome = factor(OUTCOME_SHORT[.data$outcome],
        levels = unname(OUTCOME_SHORT)
      ),
      star = dplyr::case_when(
        is.na(.data$p) ~ "", .data$p < 0.001 ~ "***",
        .data$p < 0.01 ~ "**", .data$p < 0.05 ~ "*", TRUE ~ ""
      ),
      disp = .data$value
    )

  fill_scale <- scale_fill_gradient(
    low = "grey96", high = SPEC_LEVEL_COLORS[[level]], name = "sel. freq.",
    limits = c(0, 1)
  )
  sub <- "Mean fold-selection frequency across the sparse models."

  ggplot(d, aes(.data$outcome, .data$row, fill = .data$disp)) +
    geom_tile(colour = "white", linewidth = 0.35) +
    geom_text(aes(label = .data$star),
      size = 2.1, vjust = 0.72,
      colour = "grey15", fontface = "bold"
    ) +
    facet_grid(~config, scales = "free_x", space = "free_x") +
    fill_scale +
    labs(
      x = NULL, y = NULL,
      title = sprintf(
        "%s -- %s -- across configs and outcomes",
        STAGE_LABEL[[stage]], level
      ),
      subtitle = sub
    ) +
    FIG_THEME +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4, size = 5.5),
      axis.text.y = element_text(size = 6),
      strip.text = element_text(face = "bold", size = 7),
      panel.spacing = unit(1.5, "pt"),
      plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE)
    )
}

STAGE_LABEL <- c(
  F05_classification = "classification", F06_prediction = "prediction"
)

# The F07 detail view: the strongest cells drawn from raw observations. Each
# stage shows what its claim rests on -- the feature against the phenotype, the
# held-out score against the true arm, the prediction against the truth.
pred_detail <- function(stage, level, n = 6) {
  cells <- root_cells(stage) |> filter(.data$level == !!level)
  is_class <- stage == "F05_classification"
  pcol <- if (is_class) "perm_p" else "perm_p_q2"
  cells <- cells |> slice_min(.data[[pcol]], n = n, with_ties = FALSE)
  panels <- lapply(seq_len(nrow(cells)), function(i) {
    r <- cells[i, ]
    phenotype <- if (is_class) "HR_LR" else r$outcome
    pr <- read.xlsx(
      file.path(
        sweep_root_dir(stage), level, r$config, phenotype, r$model,
        "c_data", "results.xlsx"
      ),
      "predictions"
    )
    if (!is_class) pr <- filter(pr, .data$outcome == r$outcome)
    pr$arm <- ifelse(grepl("^HR", pr$subject), "HR", "LR")
    if (is_class) {
      ggplot(pr, aes(.data$arm, .data$pred, fill = .data$arm)) +
        geom_boxplot(width = 0.5, alpha = 0.3, outlier.shape = NA) +
        geom_jitter(aes(colour = .data$arm), width = 0.12, size = 1.6) +
        scale_fill_manual(values = GROUP_COLORS, guide = "none") +
        scale_colour_manual(values = GROUP_COLORS, guide = "none") +
        labs(
          x = NULL, y = "score",
          title = sprintf(
            "%s -- %s\nAUC=%.2f p=%.3f", r$config, r$model,
            r$estimate, r$perm_p
          )
        ) +
        FIG_THEME +
        theme(plot.title = element_text(size = 6.5, face = "bold"))
    } else {
      ggplot(pr, aes(.data$y, .data$pred)) +
        geom_smooth(
          method = "lm", se = FALSE, colour = "grey30",
          linewidth = 0.5
        ) +
        geom_point(aes(colour = .data$arm), size = 1.8) +
        scale_colour_manual(values = GROUP_COLORS, guide = "none") +
        labs(
          x = r$outcome, y = "predicted",
          title = sprintf(
            "%s -- %s\nQ2=%.2f p=%.3f", r$config, r$model,
            r$q2, r$perm_p_q2
          )
        ) +
        FIG_THEME +
        theme(
          plot.title = element_text(size = 6.5, face = "bold"),
          axis.title = element_text(size = 6.5)
        )
    }
  })
  wrap_plots(panels, ncol = 3)
}

for (stage in names(STAGE_LABEL)) {
  for (level in SWEEP_LEVELS) {
    out <- file.path(REF_ROOT, stage, level)
    dir.create(out, recursive = TRUE, showWarnings = FALSE)

    hm <- stage_heatmap(stage, level)
    if (!is.null(hm)) {
      h <- if (level == "proteins") 210 else if (level == "pathways") 150 else 105
      save_panel(hm, file.path(out, "reference_heatmap"), 300, h)
    }

    detail <- pred_detail(stage, level)
    detail <- detail +
      plot_annotation(
        title = sprintf(
          "%s -- %s -- strongest cells, raw observations",
          STAGE_LABEL[[stage]], level
        ),
        subtitle = "Points by responder arm. The claim checked against the data.",
        theme = theme(
          plot.title = element_text(face = "bold", size = FIG_TITLE_SIZE),
          plot.subtitle = element_text(
            face = "italic", size = FIG_SUBTITLE_SIZE, colour = "grey30"
          )
        )
      )
    save_panel(detail, file.path(out, "reference_detail"), 230, 140)
    message(sprintf("[set] %s -- %s", stage, level))
  }
}

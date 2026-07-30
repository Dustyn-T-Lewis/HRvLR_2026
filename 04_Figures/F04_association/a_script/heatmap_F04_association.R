# Association heatmaps: one grid per feature level x method, every timepoint
# config and every outcome inside a single panel. Fill is the signed statistic,
# stars mark nominal p, and a black outline marks BH q < .05 so the corrected
# and uncorrected signals read off the same tile.

suppressWarnings(suppressMessages({
  library(here)
  source(here("functions", "sweep_composites.R"))
  source(here("functions", "sweep_drivers.R"))
  source(here("functions", "sweep_assoc_heatmap.R"))
  pacman::p_load(openxlsx, dplyr, ggplot2, patchwork)
}))

STAGE <- "F04_association"
HEATMAP_METHODS <- c("limma", "lm", "spearman")
FILL_LIMIT <- 5

assoc_heatmap <- function(level, method) {
  d <- assoc_grid(level, method)
  keep <- heatmap_features(d, level)
  keys <- driver_keys(keep, level)
  row_label <- if (level == "proteins") {
    keep
  } else {
    paste0(keys$key, ": ", keys$description)
  }
  names(row_label) <- keep

  present <- intersect(names(OUTCOME_SHORT), unique(d$outcome))
  d <- d |>
    filter(.data$feature %in% keep) |>
    mutate(
      row = factor(row_label[.data$feature], levels = rev(row_label)),
      config = factor(.data$config, levels = SWEEP_CONFIGS),
      outcome = factor(
        OUTCOME_SHORT[.data$outcome],
        levels = unname(OUTCOME_SHORT[present])
      ),
      star = case_when(
        .data$p < 0.001 ~ "***", .data$p < 0.01 ~ "**",
        .data$p < 0.05 ~ "*", TRUE ~ ""
      ),
      t_disp = pmin(pmax(.data$t, -FILL_LIMIT), FILL_LIMIT)
    )

  n_bh <- sum(d$bh < BH_CUT, na.rm = TRUE)
  ggplot(d, aes(.data$outcome, .data$row, fill = .data$t_disp)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_tile(
      data = filter(d, .data$bh < BH_CUT),
      fill = NA, colour = "black", linewidth = 0.7
    ) +
    geom_text(aes(label = .data$star),
      size = 2.4, vjust = 0.72, colour = "grey15", fontface = "bold"
    ) +
    facet_grid(~config, scales = "free_x", space = "free_x") +
    scale_fill_gradient2(
      low = "#2166AC", mid = "grey96", high = "#B2182B", midpoint = 0,
      name = "signed t", limits = c(-FILL_LIMIT, FILL_LIMIT)
    ) +
    labs(
      x = NULL, y = NULL,
      title = sprintf("%s -- %s", level, method),
      subtitle = sprintf(
        "Stars = nominal p. Black outline = BH q < .05 (%d of %d tiles).",
        n_bh, nrow(d)
      )
    ) +
    FIG_THEME +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4, size = 6),
      axis.text.y = element_text(size = 6.5),
      strip.text = element_text(face = "bold", size = 7.5),
      panel.spacing = unit(1.5, "pt"),
      plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE)
    )
}

level_composite <- function(level) {
  panels <- lapply(HEATMAP_METHODS, function(m) assoc_heatmap(level, m))
  out_dir <- file.path(sweep_root_dir(STAGE), level, "heatmap", "b_reports")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  for (i in seq_along(HEATMAP_METHODS)) {
    save_panel(
      panels[[i]],
      file.path(out_dir, paste0("heatmap_", HEATMAP_METHODS[i])), 300, 120
    )
  }

  composite <- wrap_plots(panels, ncol = 1) +
    plot_annotation(
      title = sprintf("F04 association -- %s, every config and outcome", level),
      subtitle = paste(
        "One grid per association method. Columns are the seven outcomes",
        "nested in the seven timepoint configs.\nMost tiles are pale and",
        "unstarred; outlined tiles survive within-cell BH."
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
  save_panel(composite, file.path(out_dir, "heatmap_all_methods"), 300, 330)
  message(level, " heatmaps written to ", out_dir)
}

invisible(lapply(c("modules", "pathways", "proteins"), level_composite))

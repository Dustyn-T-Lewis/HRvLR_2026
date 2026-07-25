# Association reference in the YvO F06 idiom: one feature x outcome heatmap per
# feature level, with the timepoint configs as column families. Cells carry the
# signed moderated t, stars mark nominal significance, and a bold border marks
# BH within the cell. Rows are named the way their level demands -- module ORA
# term, cleaned pathway name, gene symbol.

suppressWarnings(suppressMessages({
  library(here)
  source(here("functions", "sweep_composites.R"))
  source(here("functions", "sweep_drivers.R"))
  pacman::p_load(openxlsx, dplyr, tidyr, ggplot2, forcats, stringr, patchwork)
}))

REF_ROOT <- here("04_Figures", "shared", "reference", "b_reports")
STAGE <- "F04_association"
OUTCOMES <- c(ADAPT_OUTCOMES, "group_diff")
OUTCOME_SHORT <- c(
  comp_hypertrophy = "comp.hyp", d_fcsa_I = "fCSA I", d_fcsa_II = "fCSA II",
  d_mcsa = "mCSA", d_1rm_legpress = "1RM leg", d_1rm_ext = "1RM ext",
  group_diff = "HR-LR"
)

# Every association cell for one feature level: the moderated t and BH q per
# feature x config x outcome, taking limma as the common method so the grid is
# one model rather than a mix.
assoc_grid <- function(level, method = "limma") {
  files <- Sys.glob(file.path(
    sweep_root_dir(STAGE), level, "*", method, "c_data", "results.xlsx"
  ))
  bind_rows(lapply(files, function(f) {
    config <- basename(dirname(dirname(dirname(f))))
    sheets <- intersect(OUTCOMES, getSheetNames(f))
    bind_rows(lapply(sheets, function(s) {
      read.xlsx(f, s) |>
        transmute(
          feature = sub("@T[123]$", "", .data$feature),
          t = .data$t, p = .data$p, bh = .data$bh,
          outcome = s, config = config
        )
    }))
  })) |>
    group_by(.data$feature, .data$outcome, .data$config) |>
    slice_min(.data$p, n = 1, with_ties = FALSE) |>
    ungroup()
}

# Keep the grid legible: modules are already few, so all are shown; pathways and
# proteins are cut to the features with the strongest association anywhere.
top_features <- function(d, level, n_show) {
  if (level == "modules") {
    return(unique(d$feature))
  }
  d |>
    group_by(.data$feature) |>
    summarise(best = min(.data$p, na.rm = TRUE), .groups = "drop") |>
    slice_min(.data$best, n = n_show, with_ties = FALSE) |>
    pull(.data$feature)
}

assoc_heatmap <- function(level, n_show = 18) {
  d <- assoc_grid(level)
  keep <- top_features(d, level, n_show)
  keys <- driver_keys(keep, level)
  row_label <- if (level == "proteins") {
    keep
  } else {
    paste0(keys$key, ": ", keys$description)
  }
  names(row_label) <- keep

  d <- d |>
    filter(.data$feature %in% keep) |>
    mutate(
      row = factor(row_label[.data$feature], levels = rev(row_label)),
      config = factor(.data$config, levels = SWEEP_CONFIGS),
      outcome = factor(
        OUTCOME_SHORT[.data$outcome],
        levels = unname(OUTCOME_SHORT)
      ),
      star = case_when(
        .data$p < 0.001 ~ "***", .data$p < 0.01 ~ "**",
        .data$p < 0.05 ~ "*", TRUE ~ ""
      ),
      t_disp = pmin(pmax(.data$t, -5), 5)
    )
  hits <- filter(d, .data$bh < 0.05)

  ggplot(d, aes(.data$outcome, .data$row, fill = .data$t_disp)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_tile(
      data = hits, colour = "black", linewidth = 0.7, fill = NA
    ) +
    geom_text(aes(label = .data$star),
      size = 2.4, vjust = 0.72,
      colour = "grey15", fontface = "bold"
    ) +
    facet_grid(~config, scales = "free_x", space = "free_x") +
    scale_fill_gradient2(
      low = "#2166AC", mid = "grey96", high = "#B2182B", midpoint = 0,
      name = "signed t", limits = c(-5, 5)
    ) +
    labs(
      x = NULL, y = NULL,
      title = sprintf("%s -- association across configs and outcomes", level),
      subtitle = paste(
        "limma moderated t. Stars = nominal p; bold border = BH q<.05",
        "within cell."
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

panels <- lapply(c("modules", "pathways", "proteins"), assoc_heatmap)
composite <- wrap_plots(panels, ncol = 1, heights = c(1, 1.5, 1.5)) +
  plot_annotation(
    title = "REFERENCE -- association heatmap (YvO F06 idiom)",
    subtitle = paste(
      "Every feature level x timepoint config x outcome in one grid. Rows",
      "named per level; modules by ORA term. Read the null honestly: most",
      "cells are pale and unbordered."
    ),
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(face = "bold", size = FIG_TITLE_SIZE + 3),
      plot.subtitle = element_text(
        face = "italic", size = FIG_SUBTITLE_SIZE, colour = "grey30"
      ),
      plot.tag = element_text(face = "bold", size = FIG_TAG_SIZE)
    )
  )

dir.create(REF_ROOT, recursive = TRUE, showWarnings = FALSE)
save_panel(composite, file.path(REF_ROOT, "REFERENCE_assoc_heatmap"), 300, 330)
message("association heatmap reference written")

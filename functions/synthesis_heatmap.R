# Figure 1 of the synthesis stage: one feature x trait heatmap per level in the
# YvO F06-A idiom. Rows are named by biology and carry a colour block sized by
# the feature set it stands for. Columns are the six adaptation outcomes nested
# in the timepoint configs, each labelled with its own n. Fill is Spearman rho,
# the value is printed in the tile, stars mark nominal p and a black outline
# marks BH q < .05 applied within the cell.
#
# One method only. Spearman is fixed on design grounds -- rank-based at n =
# 13-15, no linearity or normality assumption, and the one method whose effect
# column is a single interpretable quantity across every cell. The per-method
# consequence is a separate panel, not a hidden averaging.

pacman::p_load(here, dplyr, ggplot2, patchwork, scales, stringr)
source(here("functions", "shared_style.R"))
source(here("functions", "sweep_drivers.R"))
source(here("functions", "sweep_assoc_heatmap.R"))
source(here("functions", "shared_pathway_utils.R"))

SYNTH_METHOD <- "spearman"
SYNTH_OUTCOMES <- ADAPT_OUTCOMES
SYNTH_LEVELS <- c("modules", "pathways", "proteins")
RHO_LIMIT <- 0.9

# comp_hypertrophy is the composite the HR/LR groups were cut from, so it is
# marked wherever it appears. README.md:289 bars it from any model carrying the
# HR/LR term.
CIRCULAR_OUTCOME <- "comp_hypertrophy"
SYNTH_OUTCOME_SHORT <- replace(
  OUTCOME_SHORT, CIRCULAR_OUTCOME,
  paste0(OUTCOME_SHORT[[CIRCULAR_OUTCOME]], "#")
)

# trajectory holds each feature three times, so it gets one column family per
# timepoint rather than a single tile hiding a best-of-three.
SYNTH_PANEL_CONFIGS <- c(
  "T1", "T2", "T3", "training", "acute", "total",
  "trajectory@T1", "trajectory@T2", "trajectory@T3"
)

MODULE_NAMES <- c(
  greenyellow = "Electron Transport Chain",
  magenta = "Sarcomere / Myofibril",
  black = "Aerobic Respiration",
  pink = "Cytoplasmic Translation",
  purple = "Translation (secondary)",
  yellow = "Small-Molecule Metabolism",
  brown = "Amino Acid Metabolism",
  blue = "RNA Splicing*",
  green = "",
  red = "",
  turquoise = ""
)

# Member counts come from the WGCNA atlas rather than being recounted here, so
# the label cannot drift from the network the eigengenes were built on.
module_sizes <- function() {
  atlas <- openxlsx::read.xlsx(
    here("04_Figures", "shared", "WGCNA", "c_data", "WGCNA_source_data.xlsx"),
    "module_atlas"
  )
  stats::setNames(atlas$n_proteins, atlas$module)
}

measured_proteome <- function() {
  unique(openxlsx::read.xlsx(
    here("04_Figures", "shared", "WGCNA", "c_data", "WGCNA_source_data.xlsx"),
    "module_membership"
  )$gene)
}

# Full annotated set size against the part of it this proteome actually
# measured. A strong pathway rho resting on a handful of measured members is a
# different claim from one resting on most of the set.
pathway_sizes <- function(features) {
  sets <- build_pathway_collection()
  sets <- sets[classify_database(names(sets)) %in% c("Hallmark", "GO Slim")]
  universe <- measured_proteome()
  i <- match(features, names(sets))
  data.frame(
    feature = features,
    full = lengths(sets)[i],
    measured = vapply(
      i, function(k) {
        if (is.na(k)) {
          NA_integer_
        } else {
          sum(sets[[k]] %in% universe)
        }
      },
      integer(1)
    )
  )
}

# Row identity per level: the colour block, the biology line, and the count that
# sets the block width. driver_keys already knows how each level names itself.
synthesis_rows <- function(d, level) {
  keep <- heatmap_features(d, level)
  keys <- driver_keys(keep, level)
  hits <- d |>
    filter(.data$feature %in% keep) |>
    group_by(.data$feature) |>
    summarise(
      n_hit = sum(.data$p < 0.05, na.rm = TRUE),
      n_bh = sum(.data$bh < BH_CUT, na.rm = TRUE),
      best_p = min(.data$p, na.rm = TRUE), .groups = "drop"
    )
  sizes <- if (level == "modules") module_sizes() else NULL
  out <- data.frame(
    feature = keep, key = keys$key, description = keys$description,
    fill = keys$fill
  ) |>
    left_join(hits, by = "feature")

  # Modules and pathways both carry a set size, so their block width is that
  # size; pathways stack the fraction we actually measured inside it. Proteins
  # have no size, so they get a uniform identity block sorted by evidence.
  if (level == "modules") {
    out <- out |>
      mutate(bar = sizes[.data$key], bar_inner = NA_real_) |>
      arrange(dplyr::desc(.data$bar))
  } else if (level == "pathways") {
    pw <- pathway_sizes(keep)
    out <- out |>
      mutate(
        bar = pw$full[match(.data$feature, pw$feature)],
        bar_inner = pw$measured[match(.data$feature, pw$feature)]
      ) |>
      arrange(dplyr::desc(.data$bar))
  } else {
    out <- out |>
      mutate(bar = 1, bar_inner = NA_real_) |>
      arrange(dplyr::desc(.data$n_bh), .data$best_p)
  }
  out$label <- synthesis_row_label(out, level, sizes)
  out$feature <- factor(out$feature, levels = rev(out$feature))
  out
}

synthesis_row_label <- function(rows, level, sizes) {
  if (level == "modules") {
    named <- MODULE_NAMES[rows$key]
    sprintf(
      "%s%s", str_to_title(rows$key),
      ifelse(is.na(named) | named == "", "", paste0(" - ", named))
    )
  } else if (level == "pathways") {
    sprintf("%s (%s)", gsub("\n", " ", rows$description), rows$key)
  } else {
    sprintf("%s (%s)", rows$key, sub("^module: ", "", rows$description))
  }
}

BLOCK_AXIS_TITLE <- c(
  modules = "proteins per module",
  pathways = "annotated set size (pale) vs members measured (solid)"
)

# The coloured block sits left of the grid on a reversed axis, as in the YvO
# reference. Where the feature stands for a set, the block width is that set's
# size and the name sits on the y axis so nothing clips. Proteins have no size,
# so the block is uniform and carries the name inside it.
row_block_panel <- function(rows, level) {
  sized <- level %in% c("modules", "pathways")
  p <- ggplot(rows, aes(.data$bar, .data$feature)) +
    geom_col(
      fill = rows$fill, width = 0.82, colour = "grey40", linewidth = 0.15,
      alpha = if (level == "pathways") 0.42 else 1
    ) +
    scale_x_reverse(
      expand = expansion(mult = if (sized) c(0.02, 0) else c(0.01, 0.03)),
      breaks = if (sized) scales::breaks_pretty(3) else NULL
    ) +
    labs(x = if (sized) BLOCK_AXIS_TITLE[[level]] else NULL, y = NULL) +
    FIG_THEME +
    theme(
      axis.ticks.y = element_blank(),
      axis.text.x = if (sized) element_text(size = 5.5) else element_blank(),
      axis.ticks.x = if (sized) element_line() else element_blank(),
      axis.title.x = element_text(size = 6),
      panel.grid = element_blank(),
      panel.border = element_blank()
    )
  if (level == "pathways") {
    p <- p + geom_col(
      aes(x = .data$bar_inner),
      fill = rows$fill, width = 0.82,
      colour = "grey40", linewidth = 0.15
    )
  }
  if (sized) {
    p +
      scale_y_discrete(
        labels = stats::setNames(rows$label, as.character(rows$feature))
      ) +
      theme(axis.text.y = element_text(size = 5.6, colour = "grey10"))
  } else {
    p +
      geom_text(
        aes(x = .data$bar * 0.97, label = .data$label),
        hjust = 0, size = 1.95, colour = readable_on(rows$fill)
      ) +
      theme(axis.text.y = element_blank())
  }
}

# Config strips carry their own n. n varies inside a config because the six
# outcomes are not all measured on the same subjects.
config_strip <- function(d) {
  d |>
    group_by(.data$panel_config) |>
    summarise(lo = min(.data$n), hi = max(.data$n), .groups = "drop") |>
    mutate(
      strip = ifelse(
        .data$lo == .data$hi,
        sprintf("%s\n(n=%d)", .data$panel_config, .data$lo),
        sprintf("%s\n(n=%d-%d)", .data$panel_config, .data$lo, .data$hi)
      )
    )
}

tile_grid_panel <- function(d, rows) {
  strips <- config_strip(d)
  d <- d |>
    filter(.data$feature %in% levels(rows$feature)) |>
    mutate(
      feature = factor(.data$feature, levels = levels(rows$feature)),
      outcome = factor(
        SYNTH_OUTCOME_SHORT[.data$outcome],
        levels = unname(SYNTH_OUTCOME_SHORT[SYNTH_OUTCOMES])
      ),
      strip = factor(
        strips$strip[match(.data$panel_config, strips$panel_config)],
        levels = strips$strip[
          match(SYNTH_PANEL_CONFIGS, strips$panel_config)
        ]
      ),
      star = case_when(
        .data$p < 0.001 ~ "***", .data$p < 0.01 ~ "**",
        .data$p < 0.05 ~ "*", TRUE ~ ""
      ),
      rho_disp = pmin(pmax(.data$effect, -RHO_LIMIT), RHO_LIMIT),
      shown = sprintf(
        "%.2f", ifelse(abs(.data$effect) < 0.005, 0, .data$effect)
      ),
      shown = sub("^0\\.", ".", .data$shown),
      shown = sub("^-0\\.", "-.", .data$shown)
    )
  ggplot(d, aes(.data$outcome, .data$feature, fill = .data$rho_disp)) +
    geom_tile(colour = "white", linewidth = 0.35) +
    geom_tile(
      data = filter(d, .data$bh < BH_CUT),
      fill = NA, colour = "black", linewidth = 0.65
    ) +
    geom_text(aes(label = .data$star),
      size = 2.5, vjust = -0.35, colour = "grey10", fontface = "bold"
    ) +
    geom_text(
      aes(label = .data$shown, fontface = ifelse(.data$p < 0.05, 2, 1)),
      size = 1.8, vjust = 1.15, colour = "grey10"
    ) +
    facet_grid(~strip, scales = "free_x", space = "free_x") +
    scale_fill_gradient2(
      low = "#2166AC", mid = "grey96", high = "#B2182B", midpoint = 0,
      limits = c(-RHO_LIMIT, RHO_LIMIT), name = "Spearman rho",
      breaks = c(-0.8, -0.4, 0, 0.4, 0.8)
    ) +
    labs(x = NULL, y = NULL) +
    FIG_THEME +
    theme(
      axis.text.x = element_text(
        angle = 90, hjust = 1, vjust = 0.4, size = 5.5
      ),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      strip.text = element_text(face = "bold", size = 6),
      panel.spacing = unit(1.6, "pt"),
      panel.grid = element_blank(),
      legend.position = "bottom",
      legend.key.width = unit(14, "mm"),
      legend.key.height = unit(2.6, "mm"),
      legend.title = element_text(size = 6),
      legend.text = element_text(size = 5.5)
    )
}

synthesis_subtitle <- function(d, rows, level) {
  n_pool <- length(unique(d$feature))
  n_screen <- length(unique(paste(d$config, d$outcome)))
  selected <- if (nrow(rows) < n_pool) {
    sprintf(
      "%d of %d %s shown, BH survivors first so the cap cannot drop one",
      nrow(rows), n_pool, level
    )
  } else {
    sprintf("all %d %s shown", n_pool, level)
  }
  paste0(
    sprintf(
      paste(
        "Spearman rho; %s. Screen = %d cells (%d configs x %d outcomes);",
        "trajectory is one cell drawn as its three timepoint copies. "
      ),
      selected, n_screen, length(unique(d$config)), length(unique(d$outcome))
    ),
    sprintf(
      paste(
        "%d of %d tiles reach nominal p < .05 where chance predicts %.0f,",
        "%d survive BH. "
      ),
      sum(d$p < 0.05, na.rm = TRUE), nrow(d), 0.05 * nrow(d),
      sum(d$bh < BH_CUT, na.rm = TRUE)
    ),
    "BH is applied within each cell over that cell's own feature list, ",
    "never across the screen."
  )
}

# Caveats that belong on the page rather than in a methods file: the circular
# outcome, the redundant config, the omitted contrast, and the one row whose
# name its own annotation does not support.
SYNTH_FOOTNOTES <- paste(
  "# comp_hypertrophy is the composite the HR/LR groups were defined from, so",
  "its column is circular with the grouping (README.md:289).",
  "\ntrajectory@T1/@T2/@T3 are not a fourth design: the config concatenates",
  "the three timepoint vectors into one cell, so @T1 repeats the T1 column on",
  "complete cases only (rho r = 0.96) and the three copies share a BH family.",
  "\ngroup_diff is not shown: Spearman runs the six continuous outcomes only,",
  "and the HR-vs-LR contrast is fitted in 03_DEP with subject blocking on the",
  "non-imputed matrix."
)

SYNTH_CAPTION <- c(
  modules = paste(
    "* ME_blue's top ORA term is RNA splicing, but its highest-membership",
    "proteins are pentose-phosphate enzymes (TKT, TALDO1) and leukocyte",
    "markers (LCP1, CLIC1, IFITM1). Unnamed rows returned no enriched term.",
    SYNTH_FOOTNOTES
  ),
  pathways = SYNTH_FOOTNOTES,
  proteins = paste(
    "Row label gives the protein's WGCNA module.", SYNTH_FOOTNOTES
  )
)

build_fig1_level <- function(level) {
  d <- assoc_grid(level, SYNTH_METHOD, collapse_tp = FALSE) |>
    filter(.data$outcome %in% SYNTH_OUTCOMES) |>
    mutate(
      panel_config = ifelse(
        is.na(.data$tp), .data$config, paste0(.data$config, "@", .data$tp)
      )
    )
  rows <- synthesis_rows(d, level)
  panel <- row_block_panel(rows, level) + tile_grid_panel(d, rows) +
    plot_layout(widths = c(1, 3.4)) +
    plot_annotation(
      title = sprintf(
        "%s x adaptation associations across the timepoint configs",
        str_to_sentence(level)
      ),
      subtitle = str_wrap(synthesis_subtitle(d, rows, level), 150),
      caption = str_wrap(SYNTH_CAPTION[[level]], 150),
      theme = theme(
        plot.title = element_text(face = "bold", size = FIG_TITLE_SIZE),
        plot.subtitle = element_text(
          face = "italic", size = FIG_SUBTITLE_SIZE, colour = "grey30"
        ),
        plot.caption = element_text(
          hjust = 0, size = 5.4, colour = "grey35", lineheight = 1.15
        )
      )
    )
  list(panel = panel, rows = rows, data = d)
}

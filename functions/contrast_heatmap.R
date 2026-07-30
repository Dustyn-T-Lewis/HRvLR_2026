# The HR-vs-LR contrast heatmaps, one per feature level, in the YvO F06-A idiom.
# Rows are named by biology and carry a colour block sized by the feature set
# they stand for. Columns are the nine stage 03 contrasts in three families,
# each family on its own fill scale because a difference of differences is not
# same quantity as a simple contrast. Fill is logFC, the value is printed in the
# tile, stars mark nominal p and a black outline marks BH q < .05 applied within
# the contrast.
#
# Every number here is read from 04_Features. Nothing is refitted at draw time.

pacman::p_load(
  here, dplyr, readr, ggplot2, patchwork, scales, stringr, openxlsx
)
source(here("functions", "shared_style.R"))
source(here("functions", "sweep_drivers.R"))
source(here("functions", "shared_pathway_utils.R"))

BH_CUT <- 0.05
N_SHOW <- c(modules = Inf, pathways = 18, proteins = 18)

# fgsea admits a set only when 15 of its members are measured. Nothing filters
# the singscore sets that way, so the floor is applied at row selection instead:
# HALLMARK_KRAS_SIGNALING_UP scored a BH survivor on 11 of 200 members detected.
MIN_DETECTED <- 15L

CONTRAST_FAMILIES <- list(
  list(
    id = "hrvlr", label = "HR - LR", legend = "logFC  HR - LR",
    contrasts = c(
      Baseline_HRvLR = "Baseline", Trained_HRvLR = "Trained",
      Acute_HRvLR = "Acute"
    )
  ),
  list(
    id = "interaction", label = "Interaction", legend = "logFC  interaction",
    contrasts = c(
      Training_Interaction = "Training", Acute_Interaction = "Acute"
    )
  ),
  list(
    id = "within", label = "Within-arm change", legend = "logFC  within-arm",
    contrasts = c(
      Training_HR = "Trn HR", Training_LR = "Trn LR",
      Acute_HR = "Acu HR", Acute_LR = "Acu LR"
    )
  )
)

LEVEL_BOOK <- c(
  modules = "module_contrasts.xlsx", pathways = "pathway_contrasts.xlsx",
  proteins = "protein_contrasts.xlsx"
)

contrast_grid <- function(level) {
  d <- read.xlsx(
    here("04_Features", level, "c_data", LEVEL_BOOK[[level]]), "contrasts"
  )
  d$label_key <- if (level == "proteins") d$gene else d$feature
  if (level != "pathways") {
    d$n_detected <- NA_integer_
    d$n_annotated <- NA_integer_
  }
  d |>
    filter(!is.na(.data$label_key), .data$label_key != "") |>
    select(
      "feature", "label_key", "contrast", "logFC", "p", "bh",
      "n_detected", "n_annotated"
    )
}

# Every BH survivor earns a row, then the strongest nominal features fill the
# rest. Ranking on p alone can drop a survivor, because BH is applied within a
# contrast and the contrasts hold 11 to 1,897 features.
contrast_features <- function(d, level) {
  n_show <- N_SHOW[[level]]
  if (level == "pathways") {
    d <- filter(d, .data$n_detected >= MIN_DETECTED)
  }
  if (is.infinite(n_show)) {
    return(list(keep = unique(d$feature), pool = length(unique(d$feature))))
  }
  survivors <- unique(d$feature[which(d$bh < BH_CUT)])
  ranked <- d |>
    filter(!.data$feature %in% survivors) |>
    group_by(.data$feature) |>
    summarise(best = min(.data$p, na.rm = TRUE), .groups = "drop") |>
    slice_min(.data$best,
      n = max(0, n_show - length(survivors)), with_ties = FALSE
    ) |>
    pull(.data$feature)
  list(keep = c(survivors, ranked), pool = length(unique(d$feature)))
}

module_sizes <- function() {
  atlas <- read.xlsx(
    here("04_Features", "modules", "c_data", "WGCNA_source_data.xlsx"),
    "module_atlas"
  )
  stats::setNames(atlas$n_proteins, atlas$module)
}

MODULE_NAMES <- c(
  greenyellow = "Electron Transport Chain",
  magenta = "Sarcomere / Myofibril",
  black = "Aerobic Respiration",
  pink = "Cytoplasmic Translation",
  purple = "Translation (secondary)",
  yellow = "Small-Molecule Metabolism",
  brown = "Amino Acid Metabolism",
  blue = "RNA Splicing*",
  green = "", red = "", turquoise = ""
)

# 01_Filtering scores every protein's correlation with per-sample haemoglobin.
# T3 biopsies carry twice the blood of T1 and T2, so an acute within-arm hit on
# a high-blood-correlation protein is a blood result until shown otherwise. The
# flag is the cohort's own 95th percentile, not an outside annotation.
BLOOD_QUANTILE <- 0.95

blood_tracking_genes <- function() {
  calls <- read_csv(
    here("01_Filtering", "c_data", "protein_calls.csv"),
    show_col_types = FALSE
  )
  cut <- stats::quantile(calls$blood_cor, BLOOD_QUANTILE, na.rm = TRUE)
  unique(calls$gene[which(calls$blood_cor > cut)])
}

contrast_row_label <- function(rows, level) {
  if (level == "modules") {
    named <- MODULE_NAMES[rows$key]
    sprintf(
      "%s%s", str_to_title(rows$key),
      ifelse(is.na(named) | named == "", "", paste0(" - ", named))
    )
  } else if (level == "pathways") {
    # One line only. A wrapped name overflows its bar, and the in-bar contrast
    # colour then paints the overflow white on white. The length a row can
    # actually show depends on its bar, so the trimming happens in the panel.
    gsub("\n", " ", rows$description)
  } else {
    sprintf(
      "%s (%s)%s", rows$key, sub("^module: ", "", rows$description),
      ifelse(rows$key %in% blood_tracking_genes(), "  [blood]", "")
    )
  }
}

contrast_rows <- function(d, level) {
  sel <- contrast_features(d, level)
  keep <- sel$keep
  labels <- unique(d[d$feature %in% keep, c("feature", "label_key")])
  keys <- driver_keys(labels$label_key, level)
  evidence <- d |>
    filter(.data$feature %in% keep) |>
    group_by(.data$feature) |>
    summarise(
      n_nominal = sum(.data$p < 0.05, na.rm = TRUE),
      n_bh = sum(.data$bh < BH_CUT, na.rm = TRUE),
      best_p = min(.data$p, na.rm = TRUE),
      n_detected = dplyr::first(.data$n_detected),
      n_annotated = dplyr::first(.data$n_annotated), .groups = "drop"
    )
  out <- data.frame(
    feature = labels$feature, key = keys$key,
    description = keys$description, fill = keys$fill
  ) |>
    left_join(evidence, by = "feature")

  out <- if (level == "modules") {
    sizes <- module_sizes()
    out |>
      mutate(bar = sizes[.data$key]) |>
      arrange(dplyr::desc(.data$bar))
  } else if (level == "pathways") {
    out |>
      mutate(bar = .data$n_detected) |>
      arrange(dplyr::desc(.data$bar))
  } else {
    out |>
      mutate(bar = 1) |>
      arrange(dplyr::desc(.data$n_bh), .data$best_p)
  }
  out$label <- contrast_row_label(out, level)
  out$feature <- factor(out$feature, levels = rev(out$feature))
  attr(out, "pool") <- sel$pool
  out
}

BLOCK_WIDTH <- c(modules = 1, pathways = 1.5, proteins = 1)

block_theme <- function() {
  FIG_THEME +
    theme(
      axis.ticks.y = element_blank(),
      axis.title.x = element_text(size = 6),
      panel.grid = element_blank(),
      panel.border = element_blank()
    )
}

module_block_panel <- function(rows) {
  ggplot(rows, aes(.data$bar, .data$feature)) +
    geom_col(
      fill = rows$fill, width = 0.82, colour = "grey40", linewidth = 0.15
    ) +
    scale_y_discrete(
      labels = stats::setNames(rows$label, as.character(rows$feature)),
      expand = expansion(add = 0.6)
    ) +
    scale_x_reverse(
      expand = expansion(mult = c(0.02, 0)), breaks = scales::breaks_pretty(3)
    ) +
    labs(x = "proteins per module", y = NULL) +
    block_theme() +
    theme(
      axis.text.y = element_text(size = 5.6, colour = "grey10"),
      axis.text.x = element_text(size = 5.5)
    )
}

# The Mito F03 ORA-bar idiom: bar length is how many members this proteome
# detected, the source database sits outside on the axis in its own colour, and
# the pathway name goes inside the bar, flipping outside when the bar is short.
pathway_block_panel <- function(rows) {
  # A name drawn inside its bar is painted in the fill's contrast colour, so any
  # part that overflows the bar lands white on white. Each row therefore gets a
  # character budget its own bar can hold: the widest bar fits about 36
  # characters, and a row too narrow to show 20 puts its name outside in dark
  # text at full length instead.
  span <- max(rows$bar, na.rm = TRUE)
  budget <- floor(rows$bar / (span / 36))
  inside <- budget >= 20
  shown <- vapply(
    seq_along(budget),
    function(i) {
      str_trunc(rows$label[i], if (inside[i]) budget[i] else 26L)
    },
    character(1)
  )
  ggplot(rows, aes(.data$bar, .data$feature)) +
    geom_col(
      fill = rows$fill, width = 0.82, colour = "grey40", linewidth = 0.15
    ) +
    geom_text(
      aes(x = ifelse(inside, 0, .data$bar), label = shown),
      hjust = 1.04, size = 1.9, fontface = "bold",
      colour = ifelse(inside, readable_on(rows$fill), "grey10")
    ) +
    scale_y_discrete(
      labels = stats::setNames(rows$key, as.character(rows$feature)),
      expand = expansion(add = 0.6)
    ) +
    # The slack belongs on the label side, not the heatmap side: an outside name
    # grows away from its bar towards the high-count end, and scale_x_reverse
    # maps that to the first expansion element. Putting the room on the heatmap
    # side both clipped those names and pushed the bars away from the tiles.
    scale_x_reverse(
      expand = expansion(mult = c(0.5, 0.02)),
      breaks = function(lim) {
        b <- scales::breaks_pretty(3)(c(0, max(rows$bar, na.rm = TRUE)))
        b[b >= 0]
      }
    ) +
    labs(x = "member proteins detected", y = NULL) +
    block_theme() +
    theme(
      axis.text.y = element_text(
        size = 5.6, face = "bold", colour = rev(rows$fill)
      ),
      axis.text.x = element_text(size = 5.5)
    )
}

protein_block_panel <- function(rows) {
  ggplot(rows, aes(.data$bar, .data$feature)) +
    geom_col(
      fill = rows$fill, width = 0.82, colour = "grey40", linewidth = 0.15
    ) +
    geom_text(
      aes(x = .data$bar * 0.97, label = .data$label),
      hjust = 0, size = 1.95, colour = readable_on(rows$fill)
    ) +
    scale_y_discrete(expand = expansion(add = 0.6)) +
    scale_x_reverse(expand = expansion(mult = c(0.01, 0.03))) +
    labs(x = NULL, y = NULL) +
    block_theme() +
    theme(
      axis.text.y = element_blank(), axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
}

row_block_panel <- function(rows, level) {
  switch(level,
    modules = module_block_panel(rows),
    pathways = pathway_block_panel(rows),
    proteins = protein_block_panel(rows)
  )
}

trim_zero <- function(x) sub("^(-?)0\\.", "\\1.", x)

# One panel per family so each keeps its own fill scale. Sharing a scale across
# families would put a simple contrast and a difference of differences on the
# same footing; their logFC ranges differ by a factor of two.
family_tile_panel <- function(d, rows, family, show_legend = TRUE) {
  cols <- family$contrasts
  fd <- d |>
    filter(
      .data$contrast %in% names(cols),
      .data$feature %in% levels(rows$feature)
    ) |>
    mutate(
      feature = factor(.data$feature, levels = levels(rows$feature)),
      column = factor(cols[.data$contrast], levels = unname(cols)),
      star = case_when(
        is.na(.data$p) ~ "", .data$p < 0.001 ~ "***", .data$p < 0.01 ~ "**",
        .data$p < 0.05 ~ "*", TRUE ~ ""
      ),
      shown = ifelse(
        is.na(.data$p), "-",
        trim_zero(sprintf(
          "%.2f", ifelse(abs(.data$logFC) < 0.005, 0, .data$logFC)
        ))
      )
    )
  lim <- max(abs(fd$logFC), na.rm = TRUE)
  lim <- ceiling(lim * 20) / 20
  ggplot(fd, aes(.data$column, .data$feature, fill = .data$logFC)) +
    geom_tile(colour = "grey55", linewidth = 0.3) +
    geom_tile(
      data = filter(fd, !is.na(.data$bh), .data$bh < BH_CUT),
      fill = NA, colour = "black", linewidth = 0.8
    ) +
    geom_text(
      aes(
        label = paste0(.data$shown, .data$star),
        fontface = ifelse(!is.na(.data$p) & .data$p < 0.05, 2, 1)
      ),
      size = 1.85, colour = "grey10"
    ) +
    facet_grid(~ factor(family$label, levels = family$label)) +
    scale_y_discrete(expand = expansion(add = 0.6)) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "grey96", high = "#B2182B", midpoint = 0,
      na.value = "grey85",
      limits = c(-lim, lim), breaks = c(-lim, 0, lim),
      labels = trim_zero(sprintf("%.2f", c(-lim, 0, lim))),
      name = family$legend,
      guide = if (show_legend) {
        guide_colourbar(title.position = "top", title.hjust = 0.5)
      } else {
        "none"
      }
    ) +
    labs(x = NULL, y = NULL) +
    FIG_THEME +
    theme(
      axis.text.x = element_text(
        angle = 90, hjust = 1, vjust = 0.4, size = 5.5
      ),
      axis.text.y = element_blank(), axis.ticks.y = element_blank(),
      strip.text = element_text(face = "bold", size = 6),
      panel.grid = element_blank(),
      legend.position = "bottom",
      legend.key.width = unit(11, "mm"), legend.key.height = unit(2.4, "mm"),
      legend.title = element_text(size = 5.6),
      legend.text = element_text(size = 5.2)
    )
}

contrast_subtitle <- function(d, rows, level) {
  pool <- attr(rows, "pool")
  selected <- if (nrow(rows) < pool) {
    sprintf(
      "%d of %d %s shown, BH survivors first so the cap cannot drop one",
      nrow(rows), pool, level
    )
  } else {
    sprintf("all %d %s shown", pool, level)
  }
  n_tested <- sum(!is.na(d$p))
  n_bh <- sum(d$bh < BH_CUT, na.rm = TRUE)
  paste0(
    sprintf(
      "Nine stage 03 contrasts; %s. Screen = %d estimable tests of %d cells. ",
      selected, n_tested, nrow(d)
    ),
    sprintf(
      "%d reach nominal p < .05 where chance predicts %.0f, and %s. ",
      sum(d$p < 0.05, na.rm = TRUE), 0.05 * n_tested,
      if (n_bh) sprintf("%d survive BH", n_bh) else "none survive BH"
    ),
    "BH is applied within each contrast over that contrast's own feature ",
    "list, never across the nine."
  )
}

# Caveats that belong on the page rather than in a methods file.
SHARED_FOOTNOTES <- paste(
  "# The interaction columns are the contrasts that answer the responder",
  "question: a difference at one timepoint can reflect a baseline difference,",
  "while a difference of differences cannot.",
  "\n# T3 biopsies carry roughly twice the blood of T1 and T2 (p = 0.0006).",
  "That confound cancels in the interaction but not in the two acute",
  "within-arm columns, so Acute_HR and Acute_LR are read with it in mind.",
  "\n# At 7 HR and 8 LR a nominal p < .05 needs a large effect; these columns",
  "are a declared screen, not a discovery engine.",
  "\n# A grey cell marked - is inestimable, not null: the contrast touches a",
  "Group_Time cell in which the feature was never observed, so limma returns",
  "NA and the feature is not counted in that contrast's BH denominator."
)

contrast_caption <- function(level, rows) {
  extra <- switch(level,
    modules = paste0(
      "",
      "Module detection is unsupervised and never saw the HR/LR labels, but it",
      "did see all 45 samples, and eBayes has 11 features to borrow variance",
      "across, which is close to no moderation at all. Eigengenes come from",
      "the missForest-imputed matrix, unlike the protein panel.",
      "\n* ME_blue's top ORA term is RNA splicing, but its highest-membership",
      "proteins are pentose-phosphate enzymes (TKT, TALDO1) and leukocyte",
      "markers (LCP1, CLIC1, IFITM1). Unnamed rows returned no enriched term."
    ),
    pathways = paste(
      "Rows are restricted to sets with at least 15 members detected, the same",
      "floor fgsea applies; HALLMARK_KRAS_SIGNALING_UP is excluded on that",
      "rule having scored a BH survivor on 11 of 200 annotated members. Bar",
      "length is members detected, so a short bar is a thin readout of the",
      "set.",
      "Scores",
      "come from the missForest-imputed matrix, unlike the protein panel."
    ),
    proteins = paste(
      "Read directly from stage 03: limma with duplicateCorrelation blocking",
      "on subject, fitted on the non-imputed cycloess matrix, no imputation",
      "anywhere. Row label gives the protein's WGCNA module.",
      sprintf(
        paste(
          "[blood] marks a protein whose abundance sits in the top %.0f%% of",
          "correlation with per-sample haemoglobin (01_Filtering); %d of the",
          "%d rows shown carry it, and %d are among the BH survivors."
        ),
        100 * (1 - BLOOD_QUANTILE), sum(grepl("blood", rows$label)), nrow(rows),
        sum(grepl("blood", rows$label) & rows$n_bh > 0)
      )
    )
  )
  paste(extra, SHARED_FOOTNOTES)
}

build_contrast_level <- function(level) {
  d <- contrast_grid(level)
  rows <- contrast_rows(d, level)
  shown <- filter(d, .data$feature %in% levels(rows$feature))
  tiles <- lapply(seq_along(CONTRAST_FAMILIES), function(i) {
    family_tile_panel(shown, rows, CONTRAST_FAMILIES[[i]])
  })
  n_cols <- vapply(
    CONTRAST_FAMILIES, function(f) length(f$contrasts), integer(1)
  )
  widths <- c(BLOCK_WIDTH[[level]] * 1.6, n_cols * 0.42)
  panel <- Reduce(`+`, c(list(row_block_panel(rows, level)), tiles)) +
    plot_layout(widths = widths, nrow = 1, guides = "collect") +
    plot_annotation(
      title = sprintf(
        "%s: how high responders differ from low responders",
        str_to_sentence(level)
      ),
      subtitle = str_wrap(contrast_subtitle(d, rows, level), 118),
      caption = str_wrap(contrast_caption(level, rows), 132),
      theme = theme(
        plot.title = element_text(face = "bold", size = FIG_TITLE_SIZE),
        plot.subtitle = element_text(
          face = "italic", size = FIG_SUBTITLE_SIZE, colour = "grey30"
        ),
        plot.caption = element_text(
          hjust = 0, size = 5.4, colour = "grey35", lineheight = 1.15
        ),
        legend.position = "bottom", legend.box = "horizontal"
      )
    )
  list(panel = panel, rows = rows, data = d)
}

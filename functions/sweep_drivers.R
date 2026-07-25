# Naming and bar geometry for the named-driver panels. Each feature level names
# itself differently, so the panel keeps a consistent grammar: the categorical
# key sits on the axis, the descriptive name sits inside the bar, and the fill
# carries the key. Modules are filled by their own WGCNA colour, pathways by
# source database, proteins by the level hue.

pacman::p_load(here, dplyr, ggplot2, stringr, openxlsx, scales)
source(here("functions", "shared_style.R"))
source(here("functions", "shared_pathway_utils.R"))

PROTEIN_FILL <- "#2166AC"

# WGCNA module names are literal colours, so a module bar can be filled with
# itself. Text flips to white on dark fills.
module_fill <- function(module) {
  ok <- module %in% grDevices::colors()
  ifelse(ok, module, "grey70")
}

readable_on <- function(fill) {
  rgb <- grDevices::col2rgb(fill)
  lum <- 0.299 * rgb[1, ] + 0.587 * rgb[2, ] + 0.114 * rgb[3, ]
  ifelse(lum < 140, "white", "grey10")
}

wgcna_ora <- function() {
  read.xlsx(
    here("04_Figures", "shared", "WGCNA", "c_data", "WGCNA_source_data.xlsx"),
    "module_ora"
  ) |>
    group_by(.data$module) |>
    slice_min(.data$padj, n = 1, with_ties = FALSE) |>
    ungroup() |>
    transmute(
      feature = paste0("ME_", .data$module), module = .data$module,
      term = .data$term
    )
}

# Split a feature id into the axis key and the in-bar description its level
# calls for: module colour vs its ORA term, source database vs pathway name,
# gene symbol vs nothing further.
driver_keys <- function(features, level) {
  if (level == "modules") {
    ora <- wgcna_ora()
    i <- match(features, ora$feature)
    module <- ifelse(is.na(i), sub("^ME_", "", features), ora$module[i])
    term <- ora$term[i]
    data.frame(
      key = module,
      description = ifelse(is.na(term), "no enriched term", str_trunc(term, 42)),
      fill = module_fill(module)
    )
  } else if (level == "pathways") {
    db <- classify_database(features)
    data.frame(
      key = db,
      description = str_trunc(clean_pathway_name(features, 60), 46),
      fill = unname(ifelse(db %in% names(DB_COLORS), DB_COLORS[db], "grey60"))
    )
  } else {
    # A gene symbol names itself, so the bar carries the protein's WGCNA module
    # instead: the panel then traces protein -> module without a second figure.
    mm <- read.xlsx(
      here("04_Figures", "shared", "WGCNA", "c_data", "WGCNA_source_data.xlsx"),
      "module_membership"
    )
    module <- mm$module[match(features, mm$gene)]
    data.frame(
      key = features,
      description = ifelse(is.na(module), "unassigned",
        paste("module:", module)
      ),
      fill = ifelse(is.na(module), PROTEIN_FILL, module_fill(module))
    )
  }
}

# Ranked driver bars. `score` is whatever the stage measures (selection
# frequency or -log10 p); the axis carries the key, the bar carries the name.
driver_bars <- function(d, level, xlab, title, subtitle) {
  if (!nrow(d)) {
    return(
      ggplot() +
        annotate("text", 0, 0,
          label = "no sparse signature\nin this cell",
          size = 3, colour = "grey45", fontface = "italic"
        ) +
        labs(title = title, subtitle = subtitle) +
        theme_void() +
        theme(
          plot.title = element_text(face = "bold", size = FIG_TITLE_SIZE),
          plot.subtitle = element_text(
            face = "italic", size = FIG_SUBTITLE_SIZE, colour = "grey30"
          )
        )
    )
  }
  keys <- driver_keys(d$feature, level)
  d <- d |>
    mutate(
      key = keys$key, description = keys$description, fill = keys$fill,
      text_colour = readable_on(keys$fill),
      row = stats::reorder(paste(.data$key, seq_len(dplyr::n())), .data$score)
    )

  p <- ggplot(d, aes(.data$score, .data$row)) +
    geom_col(fill = d$fill, width = 0.78, colour = "grey25", linewidth = 0.2) +
    scale_y_discrete(labels = function(x) sub(" [0-9]+$", "", x)) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.04))) +
    labs(x = xlab, y = NULL, title = title, subtitle = subtitle) +
    FIG_THEME +
    theme(
      axis.text.y = element_text(size = 7, face = "bold"),
      plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE)
    )

  if (any(nzchar(d$description))) {
    p <- p +
      geom_text(
        aes(x = 0, label = .data$description),
        hjust = -0.03, size = 2.2,
        colour = d$text_colour, fontface = "bold"
      )
  }
  p
}

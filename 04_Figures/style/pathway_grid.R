# The pathway dot-matrix (F3, S3, F4): one row per gene set within the window that selected it, one
# column per contrast or trait under its column strip, a dot per cell (fill NES in PiYG, size
# -log10 FDR) and a database square beside each row. pathway_grid() only draws; cluster_grid()
# picks the sets for F3 and S3. Expects dplyr and ggplot2 attached and figure_style.R sourced.
#
# `cells`: one row per drawn dot, with window, set_id, pathway, collection, column, nes, fdr.
# `rows`: one row per set and window, with window (a factor sets the block order), set_id and
# order (1 at the top of its block); optional block (a factor of strip labels that splits a window
# into several row strips, each in its window's colour), label (default from pathway), shade (grey
# band behind the row) and bold (label in bold).
# `columns`: column, column_window and markup (the axis label, may hold in_colour() spans), in
# display order.
pathway_collections <- c(
  Hallmark = "#332288", Reactome = "#88CCEE", KEGG_Legacy = "#DDCC77", GOBP = "#777777",
  GO_Slim = "#CC6677"
)
collection_names <- c(
  Hallmark = "Hallmark", Reactome = "Reactome", KEGG_Legacy = "KEGG", GOBP = "GO BP",
  GO_Slim = "GO Slim"
)

# The nine contrasts as columns: window strip, short label, HR and LR in their group shade.
contrast_columns <- function(contrasts) {
  tibble::tibble(
    column = contrasts, column_window = as.character(contrast_panel(contrasts)),
    markup = contrast_markup(contrasts)
  )
}

# A set's row label; `*` marks a set whose FDR < 0.05 calls in its window disagree in sign, `†` a
# training or acute set whose interaction also reaches FDR < 0.05.
set_label <- \(pathway, same_sign, interaction = FALSE) {
  paste0(
    pathway_label(pathway), if_else(same_sign, "", " *"),
    if_else(interaction %in% TRUE, " \u2020", "")
  )
}

# The candidate sets of 02_HRvLR/02_Set_Contrasts, the window as a factor in strip order (the step
# names the between-group window "HR - LR"), and each block's strip label: the step's block name,
# with the candidate count for Training and Acute patterns.
set_selection <- function(selection) {
  selection |>
    mutate(
      window = factor(if_else(window == "HR - LR", "Between groups", window), panel_families),
      block = if_else(
        is.na(pattern), as.character(window), sprintf("%s (%d)", sub(": ", "\n", block), n())
      ),
      .by = c(window, pattern)
    ) |>
    arrange(window, factor(pattern, c("shared", "HR only", "LR only"))) |>
    mutate(block = factor(block, unique(block)))
}

# One set per tree cluster at `cut` (F3, S3), the first `cap` per block by rank, in all nine
# contrasts; `windows` limits the blocks drawn. Returns the plot, the cells drawn (checked against
# set_tests) and a height in mm.
cluster_grid <- function(saved, cut, title, subtitle, cap = Inf, windows = panel_families) {
  fgsea <- filter(saved$set_tests, method == "fgsea")
  columns <- contrast_columns(unique(fgsea$contrast))
  selection <- filter(set_selection(saved$selection), window %in% windows)
  rows <- selection |>
    filter(.data[[paste0("representative_", cut)]]) |>
    arrange(rank) |>
    filter(row_number() <= cap, .by = block) |>
    left_join(distinct(fgsea, set_id, pathway), by = "set_id") |>
    mutate(order = rank, label = set_label(pathway, same_sign, interaction))
  cells <- fgsea |>
    inner_join(select(rows, window, set_id), by = "set_id", relationship = "many-to-many") |>
    transmute(window, set_id, pathway, collection, column = contrast, nes, fdr)
  panel_data <- cells |>
    inner_join(select(rows, window, block, set_id, rank, label), by = c("window", "set_id")) |>
    transmute(
      window = as.character(window), block = sub("\n", ": ", block), rank, set_id, collection,
      label, contrast = column, nes, fdr
    ) |>
    arrange(match(window, panel_families), block, rank, match(contrast, columns$column))
  source_values <- inner_join(panel_data, fgsea, by = c("set_id", "contrast"))
  clusters <- summarise(selection, n = n_distinct(.data[[paste0("cut_", cut)]]), .by = block)
  stopifnot(
    nrow(columns) == 9, all(count(panel_data, window, set_id)$n == 9),
    nrow(source_values) == nrow(panel_data),
    identical(source_values$nes.x, source_values$nes.y),
    identical(source_values$fdr.x, source_values$fdr.y),
    nrow(rows) == sum(pmin(clusters$n, cap))
  )
  list(
    panel = pathway_grid(cells, rows, columns, title = title, subtitle = subtitle),
    panel_data = panel_data, height = 3.3 * nrow(rows) + 34
  )
}

pathway_grid <- function(cells, rows, columns, title = NULL, subtitle = NULL, nes_limit = NULL,
                         fdr_cap = 12, x_angle = 0, text_size = 5.2) {
  if (!is.factor(rows$window)) rows$window <- factor(rows$window, unique(rows$window))
  pathways <- distinct(cells, set_id, pathway, collection)
  rows <- rows |>
    select(-any_of(c("pathway", "collection"))) |>
    left_join(pathways, by = "set_id") |>
    mutate(
      row = paste(window, set_id),
      block = if ("block" %in% names(rows)) block else window,
      label = if ("label" %in% names(rows)) label else pathway_label(pathway),
      shade = if ("shade" %in% names(rows)) shade else FALSE,
      bold = if ("bold" %in% names(rows)) bold else FALSE
    ) |>
    arrange(block, desc(order))
  # A narrow first column, with an empty strip, carries each row's database square.
  x_windows <- c("", unique(columns$column_window))
  x_columns <- c(" ", columns$column)
  place <- \(d) {
    mutate(d,
      row = factor(row, rows$row), block = factor(block, levels(rows$block)),
      column = factor(column, x_columns), column_window = factor(column_window, x_windows)
    )
  }
  drawn <- cells |>
    select(window, set_id, column, nes, fdr) |>
    inner_join(
      select(rows, window, block, set_id, row, collection, shade),
      by = c("window", "set_id")
    ) |>
    inner_join(select(columns, column, column_window), by = "column") |>
    mutate(log_fdr = pmin(-log10(fdr), fdr_cap)) |>
    place()
  stopifnot(nrow(drawn) == nrow(cells))
  databases <- place(mutate(rows, column = " ", column_window = ""))
  bands <- bind_rows(
    filter(drawn, shade),
    filter(databases, shade)
  )
  nes_limit <- nes_limit %||% (ceiling(max(abs(drawn$nes), na.rm = TRUE) * 10) / 10)
  # The four MSigDB databases, or GO Slim alone.
  shown <- setdiff(names(pathway_collections), "GO_Slim")
  if (all(rows$collection == "GO_Slim")) shown <- "GO_Slim"
  strips <- contrast_strips(
    x = x_windows[x_windows %in% c("", as.character(drawn$column_window))],
    y = as.character(rows$window[match(levels(droplevels(rows$block)), rows$block)])
  )
  markdown <- any(rows$bold)
  y_labels <- if (markdown) gsub("*", "\\*", rows$label, fixed = TRUE) else rows$label
  y_labels <- if_else(rows$bold, paste0("**", y_labels, "**"), y_labels)
  y_text <- if (markdown) ggtext::element_markdown else element_text
  # geom_blank() trains the axes on every row first, so the bands layer cannot reorder them.
  ggplot(drawn, aes(column, row)) +
    geom_blank() +
    geom_tile(data = bands, fill = "grey90", width = 1, height = 1) +
    # fgsea leaves a few cells without a score; they stay empty.
    geom_point(
      aes(size = log_fdr, fill = nes), filter(drawn, !is.na(nes), !is.na(fdr)),
      shape = 21, stroke = 0.2
    ) +
    # Every database keeps a key even where no row uses it, so side-by-side legends match.
    geom_point(
      aes(colour = collection), databases,
      shape = 15, size = 1.8, show.legend = c(colour = TRUE, size = FALSE, fill = FALSE)
    ) +
    scale_colour_manual(
      values = pathway_collections, limits = shown, labels = collection_names[shown],
      name = "Database"
    ) +
    scale_y_discrete(labels = rlang::set_names(y_labels, rows$row)) +
    scale_x_discrete(labels = c(` ` = "", tibble::deframe(select(columns, column, markup)))) +
    scale_fill_distiller(
      palette = "PiYG", direction = 1, limits = c(-nes_limit, nes_limit), oob = scales::squish,
      name = "NES"
    ) +
    scale_size_area(
      max_size = 2.8, limits = c(0, fdr_cap), breaks = c(2, 6, fdr_cap),
      labels = c("2", "6", paste0(fdr_cap, "+")), name = "-log10 FDR"
    ) +
    guides(
      fill = guide_colourbar(order = 1, barheight = unit(14, "mm"), barwidth = unit(2.5, "mm")),
      size = guide_legend(order = 2), colour = guide_legend(order = 3)
    ) +
    ggh4x::facet_grid2(
      block ~ column_window,
      scales = "free", space = "free", switch = "y", strip = strips
    ) +
    labs(x = NULL, y = NULL, title = title, subtitle = subtitle) +
    figure_theme +
    matrix_theme +
    theme(
      axis.text.y.left = y_text(size = text_size, hjust = 1),
      axis.text.x.bottom = ggtext::element_markdown(
        size = 5.5, angle = x_angle, hjust = if (x_angle) 1 else 0.5,
        vjust = if (x_angle == 90) 0.5 else 1
      )
    )
}

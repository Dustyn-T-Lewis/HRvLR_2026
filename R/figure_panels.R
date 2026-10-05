# Panel builders shared by the figures in 06_Figures. Each takes its data as arguments and returns
# a ggplot (or a list holding one and the numbers it plots); colours and theme come from
# figure_theme.R, which must be sourced first.

# PCA scores with a shaded 95% t-ellipse per time point and a test summary in the top-left corner.
pca_panel <- function(scores, variance, title, subtitle, stats) {
  axis <- \(i) sprintf("PC%d (%.1f%%)", i, 100 * variance[i])
  ggplot2::ggplot(scores, ggplot2::aes(PC1, PC2)) +
    ggplot2::stat_ellipse(
      ggplot2::aes(colour = timepoint, fill = timepoint),
      geom = "polygon", type = "t", level = 0.95, linewidth = 0.3, alpha = 0.12
    ) +
    ggplot2::geom_point(ggplot2::aes(colour = timepoint, shape = arm), size = 1.3) +
    ggplot2::annotate(
      "text",
      x = -Inf, y = Inf, label = stats, hjust = -0.05, vjust = 1.3,
      size = 1.9, lineheight = 0.95, colour = "grey20"
    ) +
    ggplot2::scale_colour_manual(
      values = time_colours, limits = names(time_colours), name = NULL
    ) +
    ggplot2::scale_fill_manual(values = time_colours, limits = names(time_colours), name = NULL) +
    ggplot2::scale_shape_manual(values = arm_shapes, name = NULL) +
    # Headroom above the points for the test summary.
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.32))) +
    ggplot2::labs(x = axis(1), y = axis(2), title = title, subtitle = subtitle) +
    figure_theme
}

# One window's signature rank: each arm's sets for that window (p < 0.05 or Pi < 0.05) placed in
# the other arm's ranking for the same window, drawn as limma's barcode (Wu & Smyth 2012).
# Members are ticks, the up-set above the axis and the down-set below, with a tricube moving
# average of membership along the rank. The two arms are different people, so the comparison
# shares no group mean and fry's result (02_Contrasts, fry_concordance) is a clean test.
signature_barcode <- function(proteins, fry_rows, window, set_level = "nominal",
                              span = 0.12, floor = 0.34) {
  contrasts <- paste0(window, c("_HR", "_LR"))
  threshold <- c(nominal = "p", pi = "pi_score")[[set_level]]
  level_text <- c(nominal = "p < 0.05", pi = "\u03a0 < 0.05")[[set_level]]
  rankings <- proteins |>
    dplyr::filter(contrast %in% contrasts) |>
    dplyr::mutate(rank = rank(-t, ties.method = "first"), .by = contrast) |>
    dplyr::transmute(ranked_in = sub(".*_", "", contrast), uniprot_id, rank)
  n_ranked <- dplyr::n_distinct(rankings$uniprot_id)
  lane_name <- \(chosen, ranked_in) {
    paste0(chosen, " ", tolower(window), " sets\nin ", ranked_in, "'s ranking")
  }
  lanes <- proteins |>
    dplyr::filter(contrast %in% contrasts, .data[[threshold]] < 0.05) |>
    dplyr::transmute(
      chosen = sub(".*_", "", contrast), direction = ifelse(logFC > 0, "up", "down"),
      uniprot_id, ranked_in = ifelse(chosen == "HR", "LR", "HR")
    ) |>
    dplyr::inner_join(rankings, by = c("ranked_in", "uniprot_id")) |>
    dplyr::summarise(ranks = list(rank), n = dplyr::n(), .by = c(chosen, ranked_in, direction)) |>
    dplyr::mutate(
      lane = factor(lane_name(chosen, ranked_in), lane_name(c("LR", "HR"), c("HR", "LR"))),
      side = ifelse(direction == "up", 1, -1)
    )
  ticks <- lanes |>
    dplyr::select(lane, direction, side, rank = ranks) |>
    tidyr::unnest(rank)
  worms <- lanes |>
    dplyr::mutate(worm = purrr::map(ranks, \(r) {
      smooth <- limma::tricubeMovingAverage(as.integer(seq_len(n_ranked) %in% r), span = span)
      tibble::tibble(
        rank = seq_len(n_ranked), height = floor + (1 - floor) * smooth / max(smooth)
      )
    })) |>
    dplyr::select(lane, direction, side, worm) |>
    tidyr::unnest(worm) |>
    dplyr::mutate(y = side * height, y0 = side * floor)
  stats <- fry_rows |>
    dplyr::filter(
      set_from %in% contrasts, ranked_on %in% contrasts, startsWith(set, paste0(set_level, "_"))
    ) |>
    dplyr::transmute(
      chosen = sub(".*_", "", set_from), ranked_in = sub(".*_", "", ranked_on),
      direction = sub(".*_", "", set), n_proteins, same_sign_share, fdr
    )
  labels <- dplyr::inner_join(lanes, stats, by = c("chosen", "ranked_in", "direction")) |>
    dplyr::mutate(
      text = sprintf(
        "%s (n = %d), %.0f%% same way, FDR %.2g", direction, n, 100 * same_sign_share, fdr
      ),
      x = ifelse(direction == "up", n_ranked * 0.995, n_ranked * 0.005),
      hjust = ifelse(direction == "up", 1, 0), y = side * 0.97, vjust = ifelse(side > 0, 1, 0)
    )
  # fry tested the same proteins the ticks show.
  stopifnot(nrow(labels) == 4, all(labels$n == labels$n_proteins))
  plot <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.3) +
    ggplot2::geom_segment(
      ggplot2::aes(rank, side * 0.06, xend = rank, yend = side * 0.3, colour = direction), ticks,
      linewidth = 0.2, alpha = 0.8
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(rank, ymin = y0, ymax = y, fill = direction), worms,
      alpha = 0.25
    ) +
    ggplot2::geom_line(ggplot2::aes(rank, y, colour = direction), worms, linewidth = 0.4) +
    ggplot2::geom_text(
      ggplot2::aes(x, y, label = text, hjust = hjust, vjust = vjust, colour = direction), labels,
      size = 1.9, fontface = "bold", show.legend = FALSE
    ) +
    ggplot2::facet_wrap(~lane, ncol = 1, strip.position = "left") +
    ggplot2::scale_colour_manual(
      values = c(up = direction_colours[["Up"]], down = direction_colours[["Down"]]),
      aesthetics = c("colour", "fill"), guide = "none"
    ) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.02)) +
    ggplot2::scale_y_continuous(
      limits = c(-1.05, 1.05), breaks = NULL, expand = ggplot2::expansion(0)
    ) +
    ggplot2::labs(
      x = sprintf(
        "rank in the other arm's %s contrast, up to down (n = %d)", tolower(window), n_ranked
      ),
      y = NULL, title = paste(window, "signature rank"),
      subtitle = sprintf(
        "fry: each arm's %s %s sets in the other arm's ranking", level_text, tolower(window)
      )
    ) +
    figure_theme +
    ggplot2::theme(
      strip.placement = "outside", strip.text.y.left = ggplot2::element_text(angle = 90, size = 6),
      panel.grid.major = ggplot2::element_blank(), panel.spacing.y = ggplot2::unit(2, "mm"),
      panel.border = ggplot2::element_rect(fill = NA, colour = "grey60", linewidth = 0.3)
    )
  list(
    plot = plot,
    labels = dplyr::mutate(labels, window = window, set_level = set_level, .before = 1)
  )
}

# Keeps each set unless it shares at least `cut` of its measured proteins with a set already kept;
# `ids` must be ordered strongest first.
prune_overlapping <- function(ids, members, cut = 0.5) {
  overlap <- \(a, b) length(intersect(a, b)) / min(length(a), length(b))
  purrr::reduce(ids, \(kept, id) {
    redundant <- any(purrr::map_lgl(kept, \(k) overlap(members[[id]], members[[k]]) >= cut))
    if (redundant) kept else c(kept, id)
  }, .init = character())
}

# MSigDB-style set names as short sentence-case labels with their acronyms restored.
pathway_label <- function(pathway, width = 44) {
  acronyms <- c(
    tca = "TCA", rrna = "rRNA", mrna = "mRNA", rna = "RNA", dna = "DNA", atp = "ATP",
    coa = "CoA", nmd = "NMD", myc = "MYC", rho = "RHO", gtpase = "GTPase", clec7a = "CLEC7A",
    srp = "SRP"
  )
  label <- pathway |>
    sub(pattern = "^(REACTOME|HALLMARK|KEGG|GOBP)_", replacement = "") |>
    gsub(pattern = "_", replacement = " ") |>
    stringr::str_to_sentence()
  for (word in names(acronyms)) {
    label <- gsub(paste0("\\b", word, "\\b"), acronyms[[word]], label, ignore.case = TRUE)
  }
  stringr::str_trunc(label, width)
}

# PDF (cairo, which embeds fonts and carries the Pi glyph the base device drops) and 600 dpi PNG.
save_figure <- function(plot, dir, stem, width, height) {
  ggplot2::ggsave(
    file.path(dir, paste0(stem, ".pdf")), plot,
    width = width, height = height, units = "mm", device = grDevices::cairo_pdf
  )
  ggplot2::ggsave(
    file.path(dir, paste0(stem, ".png")), plot,
    width = width, height = height, units = "mm", dpi = 600, bg = "white"
  )
}

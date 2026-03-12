# volcano_ring.R — Circular volcano-in-ring composite plot utility
# HRvLR 2026 — adapted from YvO 2025 pipeline.
# Standard Cartesian ggplot with ggforce::geom_arc_bar(); NO coord_polar().
# Ring databases: Hallmark + GO:BP + Canonical Pathways (KEGG, Reactome, WikiPathways).
# Sources pathway_utils.R for classify_database() if available.

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggforce)
  library(scales)
})

# ── Pathway utilities ─────────────────────────────────────────────────────────
# Source pathway_utils.R from YvO for classify_database() and related helpers
PW_UTILS <- file.path(dirname(rprojroot::find_rstudio_root_file()),
                      "A_YvO_2025/04_Figures_v2/shared/pathway_utils.R")
if (file.exists(PW_UTILS)) {
  source(PW_UTILS)
} else {
  # Minimal fallback classify_database if YvO not present
  classify_database <- function(pathway_names) {
    dplyr::case_when(
      grepl("^HALLMARK_", pathway_names) ~ "Hallmark",
      grepl("^REACTOME_", pathway_names) ~ "Reactome",
      grepl("^KEGG_",     pathway_names) ~ "KEGG",
      grepl("^WP_",       pathway_names) ~ "WikiPathways",
      grepl("^BIOCARTA_", pathway_names) ~ "BioCarta",
      grepl("^PID_",      pathway_names) ~ "PID",
      grepl("^GOBP_",     pathway_names) ~ "GO:BP",
      TRUE ~ "Other"
    )
  }
}

# All CP (Canonical Pathway) databases plus Hallmark and GO:BP
RING_DATABASES <- c("Hallmark", "GO:BP", "Reactome", "KEGG", "WikiPathways",
                    "BioCarta", "PID")

# ── Label cleaner ─────────────────────────────────────────────────────────────
clean_ring_label <- function(name) {
  name |>
    stringr::str_remove("^HALLMARK_") |>
    stringr::str_remove("^GOBP_") |>
    stringr::str_remove("^GOCC_") |>
    stringr::str_remove("^GOMF_") |>
    stringr::str_remove("^REACTOME_") |>
    stringr::str_remove("^KEGG_") |>
    stringr::str_remove("^WP_") |>
    stringr::str_remove("^BIOCARTA_") |>
    stringr::str_remove("^PID_") |>
    stringr::str_replace_all("_", " ") |>
    stringr::str_to_title() |>
    stringr::str_replace("Mtorc1", "mTORC1") |>
    stringr::str_replace("Myc ",   "MYC ")   |>
    stringr::str_replace("E2f ",   "E2F ")   |>
    stringr::str_replace("Dna ",   "DNA ")   |>
    stringr::str_replace("Rna ",   "RNA ")   |>
    stringr::str_replace("Tnfa ",  "TNFa ")  |>
    stringr::str_replace("Il6 ",   "IL6 ")   |>
    stringr::str_replace("Il2 ",   "IL2 ")   |>
    stringr::str_replace("Kras ",  "KRAS ")  |>
    stringr::str_replace("P53 ",   "p53 ")   |>
    stringr::str_replace("Tgf ",   "TGF ")   |>
    stringr::str_replace("Nf Kb",  "NF-kB")  |>
    stringr::str_replace("Atp ",   "ATP ")   |>
    stringr::str_replace("Ifn",    "IFN")    |>
    stringr::str_replace("Pi3k",   "PI3K")   |>
    stringr::str_replace("Akt",    "AKT")    |>
    stringr::str_replace("Mtor",   "mTOR")   |>
    stringr::str_replace("Oxidative Phosphorylation", "OXPHOS") |>
    stringr::str_replace("Mitochondrial", "Mito.")              |>
    stringr::str_replace("Organization", "Org.")                |>
    stringr::str_replace("Regulation Of", "Reg. of")           |>
    stringr::str_replace("Signaling Pathway", "Signaling")      |>
    stringr::str_replace("Biosynthetic Process", "Biosynthesis")|>
    stringr::str_replace("Catabolic Process", "Catabolism")     |>
    stringr::str_replace("Metabolic Process", "Metabolism")     |>
    stringr::str_replace("Response To", "Resp. to")             |>
    stringr::str_replace("Extracellular Matrix", "ECM")         |>
    stringr::str_wrap(width = 18)
}

# ── Ring geometry ─────────────────────────────────────────────────────────────
prepare_ring_data <- function(go_df, contrast, n_terms = 12, gap_degrees = 3,
                               start_offset = 0,
                               databases = RING_DATABASES) {
  ring <- go_df |>
    dplyr::filter(contrast == !!contrast, database %in% databases, padj < 0.05) |>
    dplyr::arrange(padj) |>
    dplyr::slice_head(n = n_terms)
  n <- nrow(ring)
  if (n == 0) {
    warning("prepare_ring_data: no significant terms for '", contrast, "'")
    return(dplyr::tibble())
  }
  arc_width_deg <- (360 - n * gap_degrees) / n
  ring |>
    dplyr::mutate(
      term_idx    = dplyr::row_number(),
      start_deg   = start_offset + (term_idx - 1) * (arc_width_deg + gap_degrees),
      end_deg     = start_deg + arc_width_deg,
      mid_deg     = (start_deg + end_deg) / 2,
      start_rad   = start_deg * pi / 180,
      end_rad     = end_deg   * pi / 180,
      mid_rad     = mid_deg   * pi / 180,
      clean_label = clean_ring_label(pathway),
      gene_list   = stringr::str_split(leadingEdge, ";")
    )
}

# ── Tick geometry ─────────────────────────────────────────────────────────────
build_tick_data <- function(ring_data, de_df, contrast,
                             tick_r0 = 4.0, tick_r1 = 5.2) {
  if (nrow(ring_data) == 0) return(dplyr::tibble())
  logfc_col <- paste0("logFC_", contrast)
  gene_lfc <- de_df |>
    dplyr::select(gene, lfc = dplyr::all_of(logfc_col)) |>
    dplyr::filter(!is.na(lfc)) |>
    dplyr::distinct(gene, .keep_all = TRUE)
  pad_rad <- 0.5 * pi / 180
  purrr::map_dfr(seq_len(nrow(ring_data)), function(i) {
    row <- ring_data[i, ]
    genes_in_arc <- intersect(row$gene_list[[1]], gene_lfc$gene)
    n_genes <- length(genes_in_arc)
    if (n_genes == 0) return(dplyr::tibble())
    arc_start <- row$start_rad + pad_rad
    arc_end   <- row$end_rad   - pad_rad
    if (arc_end <= arc_start) arc_end <- arc_start + pad_rad
    tick_angles <- seq(arc_start, arc_end, length.out = n_genes)
    matched <- gene_lfc |>
      dplyr::filter(gene %in% genes_in_arc) |>
      dplyr::slice(match(genes_in_arc, gene)) |>
      dplyr::filter(!is.na(gene))
    n_final <- nrow(matched)
    if (n_final == 0) return(dplyr::tibble())
    tick_angles <- tick_angles[seq_len(n_final)]
    dplyr::tibble(
      gene      = matched$gene,
      logFC     = matched$lfc,
      direction = ifelse(matched$lfc > 0, "Up", "Down"),
      angle_rad = tick_angles,
      x0 = tick_r0 * sin(tick_angles), y0 = tick_r0 * cos(tick_angles),
      x1 = tick_r1 * sin(tick_angles), y1 = tick_r1 * cos(tick_angles),
      term_idx  = row$term_idx,
      pathway   = row$pathway
    )
  })
}

# ── Volcano layers ────────────────────────────────────────────────────────────
build_volcano_layers <- function(de_df, contrast, volcano_radius = 3.5,
                                  fc_thresh = log2(1.5), p_thresh = 0.05,
                                  up_color   = DIR_COLORS["Up"],
                                  down_color = DIR_COLORS["Down"],
                                  ns_color   = DIR_COLORS["NS"],
                                  point_size = 0.6, point_alpha = 0.5,
                                  count_label_size = 2.8) {
  logfc_col <- paste0("logFC_", contrast)
  pval_col  <- paste0("P.Value_", contrast)
  pi_col    <- paste0("pi_score_", contrast)
  vdf <- de_df |>
    dplyr::transmute(
      gene, logFC = .data[[logfc_col]], pvalue = .data[[pval_col]],
      pi_score = .data[[pi_col]], neg_log10p = -log10(pvalue)
    ) |>
    dplyr::filter(!is.na(logFC), !is.na(pvalue), is.finite(neg_log10p)) |>
    dplyr::mutate(direction = dplyr::case_when(
      pi_score < 0.05 & logFC > 0 ~ "Up",
      pi_score < 0.05 & logFC < 0 ~ "Down",
      TRUE ~ "NS"
    ))
  n_up   <- sum(vdf$direction == "Up")
  n_down <- sum(vdf$direction == "Down")
  x_max  <- max(abs(vdf$logFC), na.rm = TRUE)
  y_max  <- max(vdf$neg_log10p, na.rm = TRUE)
  vr <- volcano_radius * 0.92
  vdf <- vdf |>
    dplyr::mutate(
      x_plot = logFC / x_max * vr,
      y_plot = (neg_log10p / y_max) * 2 * vr - vr
    )
  vdf_ns  <- vdf |> dplyr::filter(direction == "NS")
  vdf_sig <- vdf |> dplyr::filter(direction != "NS")
  layers <- list(
    ns_points  = ggplot2::geom_point(
      data = vdf_ns, ggplot2::aes(x = x_plot, y = y_plot),
      color = ns_color, size = point_size * 0.9, alpha = point_alpha * 0.6,
      inherit.aes = FALSE),
    sig_points = ggplot2::geom_point(
      data = vdf_sig, ggplot2::aes(x = x_plot, y = y_plot, color = direction),
      size = point_size * 1.3, alpha = point_alpha * 1.3, stroke = 0.3,
      inherit.aes = FALSE),
    color_scale = ggplot2::scale_color_manual(
      values = c(Up = unname(up_color), Down = unname(down_color)), guide = "none"),
    x_axis = ggplot2::annotate("segment",
      x = -vr * 0.42, xend = vr * 0.42, y = -vr, yend = -vr,
      linewidth = 0.3, linetype = "dashed", color = "grey50",
      arrow = ggplot2::arrow(ends = "both", length = ggplot2::unit(1.2, "mm"), type = "closed")),
    x_up_label = ggplot2::annotate("text",
      x = vr * 0.45, y = -vr, label = "up",
      size = count_label_size * 0.6, color = unname(up_color),
      fontface = "italic", hjust = 0),
    x_down_label = ggplot2::annotate("text",
      x = -vr * 0.45, y = -vr, label = "down",
      size = count_label_size * 0.6, color = unname(down_color),
      fontface = "italic", hjust = 1),
    x_title = ggplot2::annotate("text",
      x = 0, y = -vr - 0.35, label = "log2 FC",
      size = count_label_size * 0.6, color = "grey40", fontface = "italic"),
    y_axis = ggplot2::annotate("segment",
      x = 0, xend = 0, y = -vr, yend = vr * 0.96,
      linewidth = 0.3, linetype = "dashed", color = "grey50",
      arrow = ggplot2::arrow(ends = "last", length = ggplot2::unit(1.2, "mm"), type = "closed")),
    y_title = ggplot2::annotate("text",
      x = 0, y = vr * 1.06, label = expression(-log[10] ~ italic(p)),
      size = count_label_size * 0.6, color = "grey40"),
    n_up_box = ggplot2::annotate("label",
      x = vr * 0.5, y = vr * 0.85, label = n_up,
      size = count_label_size, fill = scales::alpha(up_color, 0.9),
      color = "black", fontface = "bold",
      label.padding = ggplot2::unit(2.5, "pt"), label.r = ggplot2::unit(2, "pt"),
      linewidth = 0.4),
    n_up_text = ggplot2::annotate("text",
      x = vr * 0.5, y = vr * 0.85, label = n_up,
      size = count_label_size, color = "white", fontface = "bold"),
    n_down_box = ggplot2::annotate("label",
      x = -vr * 0.5, y = vr * 0.85, label = n_down,
      size = count_label_size, fill = scales::alpha(down_color, 0.9),
      color = "black", fontface = "bold",
      label.padding = ggplot2::unit(2.5, "pt"), label.r = ggplot2::unit(2, "pt"),
      linewidth = 0.4),
    n_down_text = ggplot2::annotate("text",
      x = -vr * 0.5, y = vr * 0.85, label = n_down,
      size = count_label_size, color = "white", fontface = "bold")
  )
  attr(layers, "n_up")   <- n_up
  attr(layers, "n_down") <- n_down
  layers
}

# ── Ring layers ───────────────────────────────────────────────────────────────
build_ring_layers <- function(ring_data, tick_data,
                               tick_r0 = 4.0, tick_r1 = 5.2,
                               arc_r0  = 5.2, arc_r1  = 6.0,
                               up_color   = DIR_COLORS["Up"],
                               down_color = DIR_COLORS["Down"]) {
  if (nrow(ring_data) == 0) return(list())
  layers <- list()
  layers$tick_bg <- ggforce::geom_arc_bar(
    data = ring_data,
    ggplot2::aes(x0 = 0, y0 = 0, r0 = tick_r0, r = tick_r1,
                 start = start_rad, end = end_rad),
    fill = "grey93", color = "grey78", linewidth = 0.15, inherit.aes = FALSE)
  if (nrow(tick_data) > 0)
    layers$ticks <- ggplot2::geom_segment(
      data = tick_data,
      ggplot2::aes(x = x0, y = y0, xend = x1, yend = y1, color = direction),
      linewidth = 0.15, alpha = 0.7, inherit.aes = FALSE)
  layers$enrich_arcs <- ggforce::geom_arc_bar(
    data = ring_data,
    ggplot2::aes(x0 = 0, y0 = 0, r0 = arc_r0, r = arc_r1,
                 start = start_rad, end = end_rad, fill = NES),
    color = "grey40", linewidth = 0.2, inherit.aes = FALSE)
  layers$fill_scale <- ggplot2::scale_fill_gradient2(
    low = "#4393C3", mid = "white", high = "#D6604D", midpoint = 0,
    limits = c(-2, 2), oob = scales::squish, name = "NES")
  layers
}

# ── Label layer ───────────────────────────────────────────────────────────────
build_label_layer <- function(ring_data, label_r = 6.5, label_size = 2.8) {
  if (nrow(ring_data) == 0) return(list())
  label_df <- ring_data |>
    dplyr::mutate(
      arc_span   = end_rad - start_rad,
      cx         = (cos(start_rad) - cos(end_rad)) / arc_span,
      cy         = (sin(end_rad)   - sin(start_rad)) / arc_span,
      centroid_r = sqrt(cx^2 + cy^2),
      x_label    = label_r * cx / centroid_r,
      y_label    = label_r * cy / centroid_r
    )
  list(labels = ggplot2::geom_label(
    data = label_df,
    ggplot2::aes(x = x_label, y = y_label, label = clean_label),
    hjust = 0.5, vjust = 0.5, size = label_size * 0.65,
    color = "grey20", fontface = "bold", fill = "grey96",
    alpha = 0.85, lineheight = 0.8, linewidth = 0.15,
    label.padding = ggplot2::unit(1.5, "pt"),
    label.r = ggplot2::unit(0.1, "lines"),
    inherit.aes = FALSE))
}

# ── Term selection ────────────────────────────────────────────────────────────
# Picks balanced up/down terms across Hallmark + GO:BP + CP databases.
select_ring_terms <- function(go_df, contrast_name, n_each = 6,
                               databases = RING_DATABASES, min_size = 15) {
  sig <- go_df |>
    dplyr::filter(contrast == contrast_name, database %in% databases,
                  padj < 0.05, size >= min_size)

  pick_direction <- function(sig_df, n_target) {
    pool  <- sig_df |> dplyr::arrange(padj)
    hm    <- pool |> dplyr::filter(database == "Hallmark")
    bp    <- pool |> dplyr::filter(database == "GO:BP")
    cp    <- pool |> dplyr::filter(database %in% c("Reactome","KEGG","WikiPathways",
                                                     "BioCarta","PID"))
    n_hm  <- min(nrow(hm), ceiling(n_target / 3))
    n_bp  <- min(nrow(bp), ceiling(n_target / 3))
    n_cp  <- min(nrow(cp), n_target - n_hm - n_bp)
    # backfill if any database has fewer terms
    slack <- n_target - n_hm - n_bp - n_cp
    if (slack > 0) {
      extra <- head(pool |> dplyr::filter(!pathway %in% c(hm$pathway[seq_len(n_hm)],
                                                           bp$pathway[seq_len(n_bp)],
                                                           cp$pathway[seq_len(n_cp)])), slack)
      dplyr::bind_rows(
        hm |> dplyr::slice_head(n = n_hm),
        bp |> dplyr::slice_head(n = n_bp),
        cp |> dplyr::slice_head(n = n_cp),
        extra)
    } else {
      dplyr::bind_rows(
        hm |> dplyr::slice_head(n = n_hm),
        bp |> dplyr::slice_head(n = n_bp),
        cp |> dplyr::slice_head(n = n_cp))
    }
  }

  dplyr::bind_rows(
    pick_direction(sig |> dplyr::filter(NES > 0), n_each),
    pick_direction(sig |> dplyr::filter(NES < 0), n_each)
  ) |> dplyr::slice_head(n = n_each * 2)
}

# ── Ring angle helpers ────────────────────────────────────────────────────────
center_ring_angles <- function(ring, n_up) {
  n <- nrow(ring)
  if (n < 2 || n_up < 1) return(ring)
  up_mid <- (ring$start_deg[1] + ring$end_deg[min(n_up, n)]) / 2
  offset <- 90 - up_mid
  ring$start_deg <- ring$start_deg + offset
  ring$end_deg   <- ring$end_deg   + offset
  ring$mid_deg   <- ring$mid_deg   + offset
  ring$start_rad <- ring$start_deg * pi / 180
  ring$end_rad   <- ring$end_deg   * pi / 180
  ring$mid_rad   <- ring$mid_deg   * pi / 180
  ring
}

build_ring_with_gaps <- function(top_terms, contrast_name, go_df, n_each = 6,
                                  databases = RING_DATABASES) {
  real_rows  <- go_df |> dplyr::filter(contrast == contrast_name,
                                        pathway %in% top_terms$pathway)
  go_subset  <- real_rows |>
    dplyr::mutate(padj = match(pathway, top_terms$pathway) * 1e-10)
  ring <- prepare_ring_data(
    go_df = go_subset, contrast = contrast_name,
    n_terms = nrow(top_terms), gap_degrees = 3,
    start_offset = 0, databases = databases)
  n <- nrow(ring)
  if (n >= 2) {
    gap_normal <- 3; gap_split <- 8
    gaps <- rep(gap_normal, n)
    gaps[min(n_each, n)] <- gap_split
    gaps[n]              <- gap_split
    arc_width_deg <- (360 - sum(gaps)) / n
    cum_offset <- 0
    for (i in seq_len(n)) {
      if (i > 1) cum_offset <- cum_offset + arc_width_deg + gaps[i - 1]
      ring$start_deg[i] <- cum_offset
      ring$end_deg[i]   <- ring$start_deg[i] + arc_width_deg
      ring$mid_deg[i]   <- (ring$start_deg[i] + ring$end_deg[i]) / 2
      ring$start_rad[i] <- ring$start_deg[i] * pi / 180
      ring$end_rad[i]   <- ring$end_deg[i]   * pi / 180
      ring$mid_rad[i]   <- ring$mid_deg[i]   * pi / 180
    }
    ring <- center_ring_angles(ring, min(n_each, n))
  }
  ring
}

# ── Main entry point ──────────────────────────────────────────────────────────
make_volcano_ring <- function(de_df, go_df, contrast,
                               contrast_title    = NULL,
                               contrast_subtitle = NULL,
                               title_size        = 12,
                               n_terms           = 12,
                               gap_degrees       = 3,
                               start_offset      = 0,
                               databases         = RING_DATABASES,
                               volcano_radius    = 3.5,
                               tick_r0 = 4.0, tick_r1 = 5.2,
                               arc_r0  = 5.2, arc_r1  = 6.0,
                               label_r = 6.5,
                               fc_thresh = log2(1.5), p_thresh = 0.05,
                               up_color   = DIR_COLORS["Up"],
                               down_color = DIR_COLORS["Down"],
                               ns_color   = DIR_COLORS["NS"],
                               point_size = 0.6, point_alpha = 0.5,
                               label_size = 2.8, count_label_size = 2.8,
                               ring_data_override = NULL) {

  ring_data <- if (!is.null(ring_data_override)) {
    ring_data_override
  } else {
    prepare_ring_data(go_df, contrast, n_terms, gap_degrees,
                      start_offset, databases)
  }

  tick_data      <- build_tick_data(ring_data, de_df, contrast, tick_r0, tick_r1)
  volcano_layers <- build_volcano_layers(de_df, contrast, volcano_radius,
                                          fc_thresh, p_thresh,
                                          up_color, down_color, ns_color,
                                          point_size, point_alpha, count_label_size)
  ring_layers    <- build_ring_layers(ring_data, tick_data, tick_r0, tick_r1,
                                       arc_r0, arc_r1, up_color, down_color)
  label_layers   <- build_label_layer(ring_data, label_r, label_size)

  p <- ggplot2::ggplot()
  # ring background
  if (!is.null(ring_layers$tick_bg))    p <- p + ring_layers$tick_bg
  # volcano
  p <- p + volcano_layers
  # ticks + arcs
  if (!is.null(ring_layers$ticks))      p <- p + ring_layers$ticks
  if (!is.null(ring_layers$enrich_arcs))p <- p + ring_layers$enrich_arcs
  if (!is.null(ring_layers$fill_scale)) p <- p + ring_layers$fill_scale
  # labels
  if (!is.null(label_layers$labels))    p <- p + label_layers$labels

  p <- p +
    ggplot2::coord_fixed(
      xlim = c(-(label_r + 1.5), label_r + 1.5),
      ylim = c(-(label_r + 1.5), label_r + 1.8),
      clip = "off") +
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(1, 2, 2, 2, "mm"),
                   legend.position = "none")

  if (!is.null(contrast_title))
    p <- p + ggplot2::annotate("text",
      x = 0, y = label_r + 1.5, label = contrast_title,
      size = title_size / .pt, fontface = "bold", hjust = 0.5)
  if (!is.null(contrast_subtitle))
    p <- p + ggplot2::annotate("text",
      x = 0, y = label_r + 0.7, label = contrast_subtitle,
      size = (title_size - 3) / .pt, fontface = "italic",
      color = "grey30", hjust = 0.5)

  attr(p, "ring_data") <- ring_data
  attr(p, "tick_data") <- tick_data
  attr(p, "n_up")      <- attr(volcano_layers, "n_up")
  attr(p, "n_down")    <- attr(volcano_layers, "n_down")
  p
}

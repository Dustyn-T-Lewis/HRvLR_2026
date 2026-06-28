# Concordance/divergence builders shared by F04 (training) and F05 (acute).
# Each figure maps one phase's HR and LR response: agreement on the diagonal,
# responder-specific divergence off it. The divergence axis is the modelled
# interaction contrast, so the nuance rests on a test, not on the eye.

pacman::p_load(dplyr, tidyr, tibble, readr, ggplot2, ggrepel, patchwork, RRHO2, limma)

CONC_TINT <- "#FFE0E0"
DISC_TINT <- "#DCEEFF"
COMP_RED <- unname(DIR_COLORS["Up"])
COMP_BLUE <- unname(DIR_COLORS["Down"])
CLASS_COLORS <- c(
  Divergent = "#E08214", Shared = "#6A51A3",
  HR = "#2166AC", LR = "#B2182B", NS = "grey75"
)
QUAD_LEVELS <- c(
  "Concordant Up", "Concordant Down", "HR Up / LR Down", "HR Down / LR Up"
)

build_quadrant_table <- function(dep, c_hi, c_lo, c_int) {
  q <- tibble(
    gene = dep$gene,
    lfc_hi = dep[[paste0("logFC_", c_hi)]],
    lfc_lo = dep[[paste0("logFC_", c_lo)]],
    sig_hi = dep[[paste0("sig_pi_", c_hi)]] != 0,
    sig_lo = dep[[paste0("sig_pi_", c_lo)]] != 0,
    sig_int = dep[[paste0("sig_pi_", c_int)]] != 0
  ) |>
    filter(!is.na(gene), !is.na(lfc_hi), !is.na(lfc_lo)) |>
    distinct(gene, .keep_all = TRUE)

  q$quadrant <- dplyr::case_when(
    q$lfc_hi > 0 & q$lfc_lo > 0 ~ "Concordant Up",
    q$lfc_hi < 0 & q$lfc_lo < 0 ~ "Concordant Down",
    q$lfc_hi > 0 & q$lfc_lo < 0 ~ "HR Up / LR Down",
    TRUE ~ "HR Down / LR Up"
  )
  q$class <- dplyr::case_when(
    q$sig_int ~ "Divergent",
    q$sig_hi & q$sig_lo ~ "Shared",
    q$sig_hi ~ "HR",
    q$sig_lo ~ "LR",
    TRUE ~ "NS"
  )
  q
}

# Over-representation per scatter quadrant; padj_cutoff = 1 keeps non-significant
# terms so a sparse quadrant still shows its leading biology.
run_quadrant_ora <- function(quad_tbl, pw, n_show = 5) {
  universe <- quad_tbl$gene
  out <- lapply(QUAD_LEVELS, function(q) {
    genes <- quad_tbl$gene[quad_tbl$quadrant == q]
    if (length(genes) < 5) {
      return(NULL)
    }
    res <- run_ora_deduplicated(
      genes = genes, universe = universe, pathways = pw,
      jaccard_cutoff = 0.5, min_size = 15, max_size = 500, padj_cutoff = 1
    )
    if (nrow(res) == 0) {
      return(NULL)
    }
    res |>
      mutate(
        quadrant = q,
        pathway_label = clean_pathway_name(pathway, 30),
        neg_log10_padj = -log10(padj),
        significant = padj < 0.05
      ) |>
      arrange(desc(neg_log10_padj)) |>
      slice_head(n = n_show)
  })
  bind_rows(out)
}

# One quadrant's ORA as a half-bar block flanking the scatter: bars grow outward
# (left side mirrored), pathway names sit inside long bars, stars at the outer end.
make_half_bars <- function(df, fill_color, side, ylim) {
  if (is.null(df) || nrow(df) == 0) {
    return(ggplot() +
      theme_void() +
      scale_y_continuous(limits = ylim, expand = c(0, 0)))
  }
  span <- max(abs(ylim))
  n_bars <- min(nrow(df), 5L)
  y_pos <- if (ylim[1] >= 0) {
    rev(seq(0.35, span - 0.35, length.out = 5))[seq_len(n_bars)]
  } else {
    seq(-0.35, -(span - 0.35), length.out = 5)[seq_len(n_bars)]
  }
  bars <- df |>
    arrange(desc(neg_log10_padj)) |>
    slice_head(n = 5) |>
    mutate(
      y = y_pos,
      bar_fill = ifelse(significant, scales::alpha(fill_color, 0.85),
        scales::alpha(fill_color, 0.30)
      ),
      star = sig_stars(padj),
      label_inside = neg_log10_padj >= max(neg_log10_padj) * 0.5,
      label_x = ifelse(label_inside, neg_log10_padj * 0.5,
        neg_log10_padj + max(neg_log10_padj) * 0.03
      ),
      label_hjust = ifelse(label_inside, 0.5, 0),
      label_color = ifelse(label_inside,
        ifelse(significant, "white", "grey15"), "grey20"
      )
    )
  x_max <- max(bars$neg_log10_padj)
  star_mult <- if (side == "left") 0.12 else 0.04

  p <- ggplot(bars, aes(y = y)) +
    geom_rect(
      aes(xmin = 0, xmax = neg_log10_padj, ymin = y - 0.21, ymax = y + 0.21),
      fill = bars$bar_fill, color = "black", linewidth = 0.25
    ) +
    geom_text(aes(x = label_x, y = y, label = pathway_label),
      hjust = bars$label_hjust, size = 1.9, fontface = "bold",
      color = bars$label_color, lineheight = 0.85
    ) +
    geom_text(aes(x = neg_log10_padj + x_max * star_mult, label = star),
      hjust = 0, size = 2.3, fontface = "bold", color = "black"
    ) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 9) +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_blank(), axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = 6),
      axis.line.x = element_line(color = "grey60", linewidth = 0.3)
    )
  x_scale <- if (side == "left") {
    scale_x_reverse(
      limits = c(x_max * 1.18, 0), breaks = scales::pretty_breaks(3),
      expand = expansion(mult = c(0, 0))
    )
  } else {
    scale_x_continuous(
      limits = c(0, x_max * 1.18), breaks = scales::pretty_breaks(3),
      expand = expansion(mult = c(0, 0))
    )
  }
  p + x_scale +
    scale_y_continuous(limits = ylim, expand = c(0, 0)) +
    coord_cartesian(clip = "off")
}

quadrant_key <- function() {
  kd <- tibble(
    label = c("Divergent", "Shared", "HR", "LR"),
    fill = unname(CLASS_COLORS[c("Divergent", "Shared", "HR", "LR")]),
    x = c(1, 2, 2.85, 3.5)
  )
  ggplot(kd, aes(x, 0)) +
    geom_point(aes(fill = fill),
      shape = 21, size = 2.4, color = "grey45", stroke = 0.4
    ) +
    geom_text(aes(label = label), nudge_x = 0.07, hjust = 0, size = 2.2, color = "grey20") +
    scale_fill_identity() +
    scale_x_continuous(limits = c(0.5, 4.5), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-0.2, 0.2), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    theme_void()
}

# Panel A: HR-vs-LR logFC scatter (concordance map) with each quadrant's ORA as a
# half-bar block flanking its corner; interaction-significant proteins divergent.
panel_quadrant_ora <- function(quad_tbl, ora_df, labels) {
  rho <- cor(quad_tbl$lfc_hi, quad_tbl$lfc_lo, method = "spearman")
  lim <- max(abs(c(quad_tbl$lfc_hi, quad_tbl$lfc_lo))) * 1.05

  cnt <- quad_tbl |> count(quadrant, name = "total")
  sg <- quad_tbl |>
    filter(class != "NS") |>
    count(quadrant, name = "sig")
  qc <- left_join(cnt, sg, by = "quadrant") |>
    mutate(sig = tidyr::replace_na(sig, 0L))
  qlab <- function(q) {
    sprintf("%d/%d", qc$sig[qc$quadrant == q], qc$total[qc$quadrant == q])
  }
  lab_df <- quad_tbl |>
    filter(class == "Divergent") |>
    slice_max(abs(lfc_hi) + abs(lfc_lo), n = 14)

  scatter <- ggplot(quad_tbl, aes(lfc_hi, lfc_lo)) +
    annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf, fill = CONC_TINT, alpha = 0.6) +
    annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0, fill = CONC_TINT, alpha = 0.6) +
    annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0, fill = DISC_TINT, alpha = 0.6) +
    annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf, fill = DISC_TINT, alpha = 0.6) +
    geom_hline(yintercept = 0, color = "grey55", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "grey55", linewidth = 0.3) +
    geom_abline(slope = 1, linetype = "dashed", color = "black", linewidth = 0.3) +
    geom_point(data = ~ filter(.x, class == "NS"), color = "grey80", size = 0.5, alpha = 0.35) +
    geom_point(data = ~ filter(.x, class != "NS"), aes(color = class), size = 1.3, alpha = 0.9) +
    scale_color_manual(
      values = CLASS_COLORS, breaks = c("Divergent", "Shared", "HR", "LR"), name = NULL
    ) +
    ggrepel::geom_text_repel(
      data = lab_df, aes(label = gene), color = "#8A4500", size = 2.1,
      max.overlaps = 18, segment.size = 0.2, min.segment.length = 0, seed = 42
    ) +
    annotate("label",
      x = lim, y = lim, label = sprintf("Concordant Up\n%s", qlab("Concordant Up")),
      hjust = 1, vjust = 1, size = 2.4, fontface = "bold", color = COMP_RED,
      fill = alpha("white", 0.9), label.padding = unit(2, "pt"), lineheight = 0.9
    ) +
    annotate("label",
      x = -lim, y = -lim, label = sprintf("%s\nConcordant Down", qlab("Concordant Down")),
      hjust = 0, vjust = 0, size = 2.4, fontface = "bold", color = COMP_RED,
      fill = alpha("white", 0.9), label.padding = unit(2, "pt"), lineheight = 0.9
    ) +
    annotate("label",
      x = -lim, y = lim, label = sprintf("HR down / LR up\n%s", qlab("HR Down / LR Up")),
      hjust = 0, vjust = 1, size = 2.4, fontface = "bold", color = COMP_BLUE,
      fill = alpha("white", 0.9), label.padding = unit(2, "pt"), lineheight = 0.9
    ) +
    annotate("label",
      x = lim, y = -lim, label = sprintf("%s\nHR up / LR down", qlab("HR Up / LR Down")),
      hjust = 1, vjust = 0, size = 2.4, fontface = "bold", color = COMP_BLUE,
      fill = alpha("white", 0.9), label.padding = unit(2, "pt"), lineheight = 0.9
    ) +
    coord_cartesian(xlim = c(-lim, lim), ylim = c(-lim, lim), expand = FALSE) +
    labs(x = labels$x, y = labels$y) +
    FIG_THEME +
    theme(legend.position = "none", plot.margin = margin(2, 0, 0, 0, "mm"))

  nw <- make_half_bars(filter(ora_df, quadrant == "HR Down / LR Up"), COMP_BLUE, "left", c(0, lim))
  ne <- make_half_bars(filter(ora_df, quadrant == "Concordant Up"), COMP_RED, "right", c(0, lim))
  sw <- make_half_bars(filter(ora_df, quadrant == "Concordant Down"), COMP_RED, "left", c(-lim, 0))
  se <- make_half_bars(filter(ora_df, quadrant == "HR Up / LR Down"), COMP_BLUE, "right", c(-lim, 0))

  design <- c(
    patchwork::area(1, 1), patchwork::area(1, 2, 2, 2), patchwork::area(1, 3),
    patchwork::area(2, 1), patchwork::area(2, 3), patchwork::area(3, 1, 3, 3)
  )
  plot <- nw + scatter + ne + sw + se + quadrant_key() +
    plot_layout(design = design, widths = c(70, 100, 70), heights = c(85, 85, 10))
  list(plot = plot, rho = rho)
}

# Panel B: pathway NES scatter (Hallmark + GO Slim), divergent pathways from the
# cached interaction GSEA highlighted; off-diagonal distance is the magnitude gap.
build_nes_wide <- function(cache, c_hi, c_lo, c_int) {
  keep_db <- c("Hallmark", "GO Slim")
  wide <- cache |>
    filter(contrast %in% c(c_hi, c_lo), database %in% keep_db, is.finite(NES)) |>
    select(pathway, database, contrast, NES, padj, size) |>
    pivot_wider(
      id_cols = c(pathway, database),
      names_from = contrast, values_from = c(NES, padj, size)
    ) |>
    filter(
      !is.na(.data[[paste0("NES_", c_hi)]]),
      !is.na(.data[[paste0("NES_", c_lo)]])
    )
  int <- cache |>
    filter(contrast == c_int, database %in% keep_db) |>
    select(pathway, NES_int = NES, padj_int = padj)
  left_join(wide, int, by = "pathway")
}

panel_nes_scatter <- function(nes_wide, c_hi, c_lo, labels) {
  x <- nes_wide[[paste0("NES_", c_hi)]]
  y <- nes_wide[[paste0("NES_", c_lo)]]
  ct <- suppressWarnings(cor.test(x, y, method = "spearman"))
  lim <- max(abs(c(x, y))) * 1.05
  d <- nes_wide |>
    mutate(
      nx = x, ny = y,
      class = if_else(!is.na(padj_int) & padj_int < 0.05, "Divergent", "Concordant")
    )
  lab_d <- d |>
    filter(class == "Divergent") |>
    slice_max(abs(nx - ny), n = 14)
  tints <- tibble(
    xmin = c(0, -lim, 0, -lim), xmax = c(lim, 0, lim, 0),
    ymin = c(0, -lim, -lim, 0), ymax = c(lim, 0, 0, lim),
    fill = c(CONC_TINT, CONC_TINT, DISC_TINT, DISC_TINT)
  )

  p <- ggplot(d, aes(nx, ny)) +
    geom_rect(
      data = tints,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
      inherit.aes = FALSE, alpha = 0.5
    ) +
    scale_fill_identity() +
    geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "grey60", linewidth = 0.3) +
    geom_abline(slope = 1, linetype = "dashed", color = "grey45", linewidth = 0.3) +
    geom_point(aes(shape = database, color = class), size = 2, alpha = 0.85) +
    scale_shape_manual(values = c(Hallmark = 17, "GO Slim" = 16), name = NULL) +
    scale_color_manual(
      values = c(Divergent = "#E08214", Concordant = "grey55"), name = NULL
    ) +
    ggrepel::geom_text_repel(
      data = lab_d, aes(label = clean_pathway_name(pathway, 30)),
      color = "#8A4500", size = 2.1, max.overlaps = 18, segment.size = 0.2
    ) +
    coord_fixed(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
    labs(
      x = sprintf("%s NES", labels$x_short), y = sprintf("%s NES", labels$y_short),
      subtitle = sprintf(
        "Spearman rho = %.2f; %d pathways (Hallmark + GO Slim)",
        unname(ct$estimate), nrow(d)
      )
    ) +
    FIG_THEME

  list(plot = p, rho = unname(ct$estimate), data = d)
}

# GSEA-style running enrichment of a member set along a ranked statistic.
running_es <- function(t_vals, in_set) {
  n <- length(t_vals)
  n_h <- sum(in_set)
  if (n_h == 0) {
    return(rep(0, n))
  }
  hit_w <- ifelse(in_set, abs(t_vals), 0)
  cumsum(ifelse(in_set, hit_w / sum(hit_w), -1 / (n - n_h)))
}

# Supp: fry rotation test of the HR-significant up/down DEP sets along the LR
# moderated-t ranking. fry carries the within-subject block + duplicateCorrelation
# (camera cannot), so it is valid for this repeated-measures design. Returns the
# fry result, the ranked list with running ES, and leading-edge ORA per direction.
run_fry_concordance <- function(da, dep, pw, c_hi, c_lo, lo_levels) {
  set.seed(42)
  mat <- as.matrix(da$data)
  meta <- da$metadata
  meta$Group_Time <- factor(meta$Group_Time)
  design <- model.matrix(~ 0 + Group_Time, data = meta)
  colnames(design) <- sub("^Group_Time", "", colnames(design))
  block <- meta$Subject_ID
  corfit <- duplicateCorrelation(mat, design, block = block)
  cm <- makeContrasts(contrasts = paste0(lo_levels[2], " - ", lo_levels[1]), levels = design)

  imp_ids <- rownames(mat)
  sig <- dep[dep[[paste0("sig_pi_", c_hi)]] != 0 & dep$uniprot_id %in% imp_ids, ]
  lfc_hi <- sig[[paste0("logFC_", c_hi)]]
  up_ids <- sig$uniprot_id[lfc_hi > 0]
  dn_ids <- sig$uniprot_id[lfc_hi < 0]
  idx <- list(Up = match(up_ids, imp_ids), Down = match(dn_ids, imp_ids))
  idx <- idx[vapply(idx, length, integer(1)) >= 3]
  res <- fry(mat,
    index = idx, design = design, contrast = cm[, 1],
    block = block, correlation = corfit$consensus.correlation
  )
  res <- tibble::rownames_to_column(as.data.frame(res), "set")
  res$consistent <- res$Direction == res$set
  res$n <- vapply(idx[res$set], length, integer(1))

  rk <- tibble(
    gene = dep$gene, uid = dep$uniprot_id, t_lo = dep[[paste0("t_", c_lo)]]
  ) |>
    filter(!is.na(t_lo), uid %in% imp_ids) |>
    arrange(desc(t_lo)) |>
    mutate(rank = row_number(), in_up = uid %in% up_ids, in_down = uid %in% dn_ids)
  rk$es_up <- running_es(rk$t_lo, rk$in_up)
  rk$es_down <- running_es(rk$t_lo, rk$in_down)

  ora_le <- function(genes) {
    if (length(genes) < 5) {
      return(NULL)
    }
    r <- run_ora_deduplicated(genes, rk$gene, pw,
      jaccard_cutoff = 0.5, min_size = 15, max_size = 500, padj_cutoff = 1
    )
    if (nrow(r) == 0) {
      return(NULL)
    }
    r |>
      mutate(
        pathway_label = clean_pathway_name(pathway, 28),
        neg_log10_padj = -log10(padj), significant = padj < 0.05
      ) |>
      arrange(desc(neg_log10_padj)) |>
      slice_head(n = 5)
  }
  le_up <- rk$gene[rk$in_up & rk$rank <= rk$rank[which.max(rk$es_up)]]
  le_dn <- rk$gene[rk$in_down & rk$rank >= rk$rank[which.min(rk$es_down)]]

  list(
    results = res, rank = rk, ora_up = ora_le(le_up), ora_down = ora_le(le_dn),
    n_all = nrow(rk)
  )
}

ora_bars <- function(ora_df, fill, n = 5) {
  if (is.null(ora_df) || nrow(ora_df) == 0) {
    return(ggplot() +
      theme_void())
  }
  d <- ora_df |>
    arrange(desc(neg_log10_padj)) |>
    slice_head(n = n) |>
    mutate(
      term = reorder(pathway_label, neg_log10_padj),
      bar_fill = ifelse(significant, scales::alpha(fill, 0.85), scales::alpha(fill, 0.3))
    )
  ggplot(d, aes(neg_log10_padj, term)) +
    geom_col(fill = d$bar_fill, color = "black", linewidth = 0.25, width = 0.7) +
    geom_text(aes(label = sig_stars(padj)), hjust = -0.2, size = 2.3) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(x = expression(-log[10](p[adj])), y = NULL) +
    FIG_THEME +
    theme(axis.text.y = element_text(size = 6.5))
}

# Assemble the running-ES curves, barcode rugs, ranking metric, and leading-edge
# ORA into the YvO fry layout (ES + barcode per direction stacked, ORA to the side).
panel_fry <- function(fry_out, labels) {
  rk <- fry_out$rank
  res <- fry_out$results
  n_all <- fry_out$n_all
  fr <- function(s) res[res$set == s, , drop = FALSE]

  es_plot <- function(es_col, color, title, set) {
    row <- fr(set)
    is_sig <- nrow(row) && !is.na(row$PValue) && row$PValue < 0.05
    col <- if (is_sig) color else scales::alpha(color, 0.4)
    lab <- if (nrow(row)) {
      sprintf("fry %s, %s (n = %d)", set, fmt_p(row$PValue), row$n)
    } else {
      set
    }
    ggplot(rk, aes(rank, .data[[es_col]])) +
      geom_area(fill = scales::alpha(col, 0.15), color = NA) +
      geom_line(color = col, linewidth = 0.6) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.3) +
      annotate("text",
        x = n_all, y = Inf, label = lab, hjust = 1.05, vjust = 1.6,
        size = 2.4, fontface = "bold", color = "grey20"
      ) +
      labs(title = title, x = NULL, y = "ES") +
      scale_x_continuous(limits = c(1, n_all), expand = c(0.005, 0)) +
      FIG_THEME +
      theme(
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.title = element_text(size = 9, face = "bold")
      )
  }
  bc_plot <- function(in_col, color, set) {
    row <- fr(set)
    is_sig <- nrow(row) && !is.na(row$PValue) && row$PValue < 0.05
    col <- if (is_sig) color else scales::alpha(color, 0.4)
    ggplot(rk[rk[[in_col]], ], aes(x = rank, xend = rank, y = 0, yend = 1)) +
      geom_segment(color = col, linewidth = 0.3, alpha = 0.7) +
      scale_x_continuous(limits = c(1, n_all), expand = c(0.005, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      theme_void() +
      theme(
        panel.background = element_rect(fill = "grey97", color = NA),
        plot.margin = margin(0, 1, 0, 2, "mm")
      )
  }
  rank_curve <- ggplot(rk, aes(rank, t_lo)) +
    geom_area(fill = scales::alpha(COMP_BLUE, 0.2), color = COMP_BLUE, linewidth = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
    labs(x = sprintf("Rank (%s moderated t, n = %d)", labels$lo, n_all), y = NULL) +
    scale_x_continuous(limits = c(1, n_all), expand = c(0.005, 0)) +
    FIG_THEME

  design <- c(
    patchwork::area(1, 1), patchwork::area(2, 1), patchwork::area(3, 1),
    patchwork::area(4, 1), patchwork::area(5, 1),
    patchwork::area(1, 2, 2, 2), patchwork::area(3, 2, 4, 2)
  )
  es_plot("es_up", COMP_RED, sprintf("HR-up DEPs along %s rank", labels$lo), "Up") +
    bc_plot("in_up", COMP_RED, "Up") +
    es_plot("es_down", COMP_BLUE, sprintf("HR-down DEPs along %s rank", labels$lo), "Down") +
    bc_plot("in_down", COMP_BLUE, "Down") +
    rank_curve +
    ora_bars(fry_out$ora_up, COMP_RED) +
    ora_bars(fry_out$ora_down, COMP_BLUE) +
    plot_layout(
      design = design, heights = c(2, 0.25, 2, 0.25, 0.7), widths = c(1.7, 1.5)
    )
}

# Supp: RRHO2 threshold-free overlap of the two moderated-t rankings.
panel_rrho2 <- function(dep, c_hi, c_lo, labels) {
  rr <- tibble(
    gene = dep$gene,
    t1 = dep[[paste0("t_", c_hi)]],
    t2 = dep[[paste0("t_", c_lo)]]
  ) |>
    filter(!is.na(gene), !is.na(t1), !is.na(t2)) |>
    distinct(gene, .keep_all = TRUE)

  obj <- RRHO2_initialize(
    data.frame(gene = rr$gene, score = rr$t1),
    data.frame(gene = rr$gene, score = rr$t2),
    labels = c(labels$hi, labels$lo), log10.ind = TRUE,
    multipleTesting = "none", boundary = 0.02, method = "hyper", stepsize = 20
  )
  hm <- obj$hypermat
  hm_df <- expand.grid(i = seq_len(nrow(hm)), j = seq_len(ncol(hm)))
  hm_df$value <- as.vector(hm)

  plot <- ggplot(hm_df, aes(i, j, fill = value)) +
    geom_raster() +
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
      name = expression(signed ~ -log[10] ~ P)
    ) +
    annotate("text",
      x = nrow(hm) * 0.85, y = ncol(hm) * 0.9, label = "Conc. up",
      size = 2.4, fontface = "bold"
    ) +
    annotate("text",
      x = nrow(hm) * 0.15, y = ncol(hm) * 0.1, label = "Conc. down",
      size = 2.4, fontface = "bold"
    ) +
    coord_fixed() +
    labs(
      x = sprintf("%s rank", labels$hi), y = sprintf("%s rank", labels$lo),
      subtitle = sprintf("%d shared proteins", nrow(rr))
    ) +
    FIG_THEME +
    theme(
      axis.text = element_blank(), axis.ticks = element_blank(),
      panel.grid = element_blank()
    )
  list(plot = plot, n_shared = nrow(rr))
}

# Supp: how the quadrant split holds up as the significance threshold changes.
panel_threshold_sens <- function(dep, c_hi, c_lo) {
  base <- tibble(
    lfc_hi = dep[[paste0("logFC_", c_hi)]],
    lfc_lo = dep[[paste0("logFC_", c_lo)]],
    pi_hi = dep[[paste0("sig_pi_", c_hi)]],
    pi_lo = dep[[paste0("sig_pi_", c_lo)]],
    fdr_hi = dep[[paste0("adj.P.Val_", c_hi)]],
    fdr_lo = dep[[paste0("adj.P.Val_", c_lo)]],
    p_hi = dep[[paste0("P.Value_", c_hi)]],
    p_lo = dep[[paste0("P.Value_", c_lo)]]
  ) |>
    filter(!is.na(lfc_hi), !is.na(lfc_lo))
  base$quad <- dplyr::case_when(
    base$lfc_hi > 0 & base$lfc_lo > 0 ~ "Concordant Up",
    base$lfc_hi < 0 & base$lfc_lo < 0 ~ "Concordant Down",
    TRUE ~ "Discordant"
  )

  thr <- list(
    "pi (signed)" = base$pi_hi != 0 | base$pi_lo != 0,
    "FDR < 0.05" = base$fdr_hi < 0.05 | base$fdr_lo < 0.05,
    "FDR < 0.10" = base$fdr_hi < 0.10 | base$fdr_lo < 0.10,
    "Nom. p < 0.05" = base$p_hi < 0.05 | base$p_lo < 0.05
  )
  counts <- bind_rows(lapply(names(thr), function(nm) {
    base[thr[[nm]] & !is.na(thr[[nm]]), ] |>
      count(quad, name = "n") |>
      mutate(threshold = nm)
  }))
  counts$threshold <- factor(counts$threshold, levels = names(thr))

  plot <- ggplot(counts, aes(threshold, n, fill = quad)) +
    geom_col(position = "stack") +
    geom_text(aes(label = n),
      position = position_stack(vjust = 0.5), size = 2.6, color = "white"
    ) +
    scale_fill_manual(
      values = c(
        "Concordant Up" = "#E57373", "Concordant Down" = "#64B5F6",
        "Discordant" = "#FFB74D"
      ), name = "Quadrant"
    ) +
    labs(x = NULL, y = "Proteins selected") +
    FIG_THEME +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  list(plot = plot, data = counts)
}

# Supp: bootstrap CI on the HR-vs-LR logFC Spearman concordance.
panel_rho_bootstrap <- function(dep, c_hi, c_lo, reps = 1000) {
  set.seed(42)
  d <- tibble(
    x = dep[[paste0("logFC_", c_hi)]],
    y = dep[[paste0("logFC_", c_lo)]]
  ) |>
    filter(!is.na(x), !is.na(y))
  obs <- cor(d$x, d$y, method = "spearman")
  vals <- vapply(seq_len(reps), function(i) {
    idx <- sample(nrow(d), replace = TRUE)
    cor(d$x[idx], d$y[idx], method = "spearman")
  }, numeric(1))
  ci <- quantile(vals, c(0.025, 0.975))

  plot <- ggplot(tibble(rho = vals), aes(rho)) +
    geom_histogram(bins = 40, fill = "grey75", color = "white") +
    annotate("rect",
      xmin = ci[1], xmax = ci[2], ymin = 0, ymax = Inf, alpha = 0.15, fill = "#2166AC"
    ) +
    geom_vline(xintercept = obs, linetype = "dashed", color = "#B2182B", linewidth = 0.6) +
    labs(
      x = "Spearman rho (logFC HR vs LR)", y = "Bootstrap reps",
      subtitle = sprintf("rho = %.2f [%.2f, %.2f], %d reps", obs, ci[1], ci[2], reps)
    ) +
    FIG_THEME
  list(plot = plot, data = tibble(
    stat = "spearman_rho", observed = obs,
    ci_lo = ci[1], ci_hi = ci[2], reps = reps
  ))
}

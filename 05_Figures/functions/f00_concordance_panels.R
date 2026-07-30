# Concordance builders shared by the training and acute concordance leaves. Each figure maps
# one phase's HR and LR response as a 5-panel composite matching the YvO engine:
#   A quadrant ORA scatter, B pattern heatmap + Sankey, C fry rotation test,
#   D pathway NES scatter, E RRHO2. Both figures are the concordance flavour
#   (ref_slope = +1, diagonal = concordant); divergence is the modelled
#   interaction contrast, not slope. Only the per-figure cfg differs.

pacman::p_load(
  dplyr, tidyr, tibble, readr, stringr, ggplot2, ggrepel, patchwork,
  scales, cowplot, RRHO2, limma, fgsea
)

# Palette

CONC_TINT <- "#FFE0E0" # concordant diagonal quadrants
DISC_TINT <- "#DCEEFF" # discordant off-diagonal quadrants
COMP_RED <- unname(DIR_COLORS["Up"])
COMP_BLUE <- unname(DIR_COLORS["Down"])

# Panel-A significance classes: interaction-significant proteins are Divergent.
SIG_COLORS_CONC <- c(
  "Divergent"   = "#7B5EA7",
  "Sig Both"    = "#2E7D32",
  "Sig HR only" = unname(GROUP_COLORS["HR"]),
  "Sig LR only" = unname(GROUP_COLORS["LR"]),
  "NS"          = "grey70"
)
SIG_LABEL_FILL_CONC <- scales::alpha(SIG_COLORS_CONC, 0.78)
names(SIG_LABEL_FILL_CONC) <- names(SIG_COLORS_CONC)

# Panel-B pattern classes
QUAD_ORDER_B <- c("Concordant Up", "Concordant Down", "Discordant")
QUAD_COLORS_B <- c(
  "Concordant Up" = unname(GROUP_COLORS["LR"]),
  "Concordant Down" = unname(GROUP_COLORS["HR"]),
  "Discordant" = "#1B7837"
)
QUAD_BG_B <- c(
  "Concordant Up" = "#F4D9D2", "Concordant Down" = "#D5DEEF",
  "Discordant" = "#C8E0CD"
)
ENDPOINT_COLORS_B <- c(
  "Concordant Up" = "#67001F", "Concordant Down" = "#053061",
  "Discordant" = "#00441B"
)
SIG_COLORS_B <- c(
  "Both" = "#2E7D32", "HR" = unname(GROUP_COLORS["HR"]),
  "LR" = unname(GROUP_COLORS["LR"]), "Inter." = "#7B5EA7", "NS" = "grey70"
)

# GO Slim consolidation (ported from YvO F04-F06_go_slim_terms.R)

BP_SLIM <- c(
  "GO:0000278", "GO:0000910", "GO:0002181", "GO:0002376", "GO:0003012",
  "GO:0003013", "GO:0003014", "GO:0003016", "GO:0005975", "GO:0006091",
  "GO:0006260", "GO:0006281", "GO:0006310", "GO:0006325", "GO:0006351",
  "GO:0006355", "GO:0006399", "GO:0006457", "GO:0006520", "GO:0006629",
  "GO:0006766", "GO:0006886", "GO:0006913", "GO:0006914", "GO:0006954",
  "GO:0007005", "GO:0007010", "GO:0007018", "GO:0007031", "GO:0007059",
  "GO:0007126", "GO:0007155", "GO:0007163", "GO:0007586", "GO:0009100",
  "GO:0012501", "GO:0016071", "GO:0016192", "GO:0023052", "GO:0030154",
  "GO:0030163", "GO:0030198", "GO:0032200", "GO:0034330", "GO:0042060",
  "GO:0042180", "GO:0042254", "GO:0044782", "GO:0048856", "GO:0048870",
  "GO:0050877", "GO:0051604", "GO:0055085", "GO:0055086", "GO:0061024",
  "GO:0065003", "GO:0071941", "GO:0072659", "GO:0098542", "GO:0098754",
  "GO:0140014", "GO:1901135"
)

SLIM_CONSOLIDATED <- c(
  "GO:0003012" = "Muscle & Contractile",
  "GO:0003013" = "Circulatory System",
  "GO:0030198" = "ECM & Adhesion", "GO:0007155" = "ECM & Adhesion",
  "GO:0034330" = "ECM & Adhesion", "GO:0042060" = "ECM & Adhesion",
  "GO:0007010" = "Cytoskeleton & Motility", "GO:0048870" = "Cytoskeleton & Motility",
  "GO:0007018" = "Cytoskeleton & Motility", "GO:0007163" = "Cytoskeleton & Motility",
  "GO:0044782" = "Cytoskeleton & Motility",
  "GO:0002376" = "Immune & Inflammation", "GO:0006954" = "Immune & Inflammation",
  "GO:0098542" = "Immune & Inflammation",
  "GO:0006629" = "Lipid Metabolism", "GO:0042180" = "Lipid Metabolism",
  "GO:0005975" = "Carbohydrate & Energy Metabolism",
  "GO:0006091" = "Carbohydrate & Energy Metabolism",
  "GO:1901135" = "Carbohydrate & Energy Metabolism",
  "GO:0006520" = "Amino Acid & Cofactor Metabolism",
  "GO:0055086" = "Amino Acid & Cofactor Metabolism",
  "GO:0006766" = "Amino Acid & Cofactor Metabolism",
  "GO:0071941" = "Amino Acid & Cofactor Metabolism",
  "GO:0098754" = "Amino Acid & Cofactor Metabolism",
  "GO:0007586" = "Amino Acid & Cofactor Metabolism",
  "GO:0007005" = "Mitochondria & Energy", "GO:0007031" = "Mitochondria & Energy",
  "GO:0006457" = "Protein Homeostasis", "GO:0030163" = "Protein Homeostasis",
  "GO:0006914" = "Protein Homeostasis", "GO:0051604" = "Protein Homeostasis",
  "GO:0065003" = "Protein Homeostasis", "GO:0009100" = "Protein Homeostasis",
  "GO:0055085" = "Transport", "GO:0016192" = "Transport",
  "GO:0006886" = "Transport", "GO:0006913" = "Transport",
  "GO:0072659" = "Transport", "GO:0061024" = "Transport",
  "GO:0002181" = "Translation & Ribosome", "GO:0042254" = "Translation & Ribosome",
  "GO:0006399" = "Translation & Ribosome",
  "GO:0006351" = "Transcription & Chromatin", "GO:0006355" = "Transcription & Chromatin",
  "GO:0016071" = "Transcription & Chromatin", "GO:0006325" = "Transcription & Chromatin",
  "GO:0006281" = "DNA & Cell Cycle", "GO:0006260" = "DNA & Cell Cycle",
  "GO:0006310" = "DNA & Cell Cycle", "GO:0032200" = "DNA & Cell Cycle",
  "GO:0000278" = "DNA & Cell Cycle", "GO:0140014" = "DNA & Cell Cycle",
  "GO:0007059" = "DNA & Cell Cycle", "GO:0000910" = "DNA & Cell Cycle",
  "GO:0007126" = "DNA & Cell Cycle",
  "GO:0048856" = "Development", "GO:0030154" = "Development",
  "GO:0012501" = "Development", "GO:0003014" = "Development",
  "GO:0003016" = "Development"
)

CONSOLIDATED_PATHWAY_ORDER <- c(
  "Muscle & Contractile", "Cytoskeleton & Motility", "ECM & Adhesion",
  "Lipid Metabolism", "Carbohydrate & Energy Metabolism",
  "Amino Acid & Cofactor Metabolism",
  "Mitochondria & Energy", "Protein Homeostasis",
  "Transport", "Translation & Ribosome", "Transcription & Chromatin",
  "Immune & Inflammation", "DNA & Cell Cycle", "Circulatory System",
  "Development", "Other"
)

# Map foreground genes to a single consolidated GO-Slim BP category
# (ported from YvO assign_go_slim_consolidated). Each gene is assigned the
# rarest slim term it carries (Development de-prioritised as a catch-all).
assign_go_slim_consolidated <- function(fg_genes, all_genes, min_cat_size = 2) {
  requireNamespace("GO.db", quietly = TRUE)
  requireNamespace("org.Hs.eg.db", quietly = TRUE)
  requireNamespace("AnnotationDbi", quietly = TRUE)

  suppressMessages({
    all_entrez <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
      keys = all_genes, keytype = "SYMBOL", column = "ENTREZID",
      multiVals = "first"
    )
    all_go <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
      keys = as.character(stats::na.omit(all_entrez)),
      keytype = "ENTREZID", columns = c("SYMBOL", "GO", "ONTOLOGY")
    )
  })

  all_bp <- all_go |>
    filter(ONTOLOGY == "BP", !is.na(GO)) |>
    distinct(SYMBOL, GO)

  ancestors <- as.list(GO.db::GOBPANCESTOR)
  all_go_ids <- unique(all_bp$GO)

  go_to_slim <- setNames(
    lapply(all_go_ids, function(go_id) {
      hits <- character(0)
      if (go_id %in% BP_SLIM) hits <- go_id
      anc <- ancestors[[go_id]]
      if (!is.null(anc)) hits <- c(hits, intersect(anc, BP_SLIM))
      unique(hits)
    }),
    all_go_ids
  )

  all_gene_slim <- all_bp |>
    mutate(slim_list = go_to_slim[GO]) |>
    unnest(slim_list) |>
    select(SYMBOL, slim = slim_list) |>
    distinct()

  fg_gene_slim <- all_gene_slim |> filter(SYMBOL %in% fg_genes)
  fg_gene_consolidated <- fg_gene_slim |>
    mutate(consolidated = SLIM_CONSOLIDATED[slim]) |>
    filter(!is.na(consolidated))

  fg_term_counts <- fg_gene_slim |> count(slim, name = "n_fg")

  best_consolidated <- fg_gene_consolidated |>
    left_join(fg_term_counts, by = "slim") |>
    mutate(priority = ifelse(consolidated == "Development", 2, 1)) |>
    arrange(priority, n_fg) |>
    group_by(SYMBOL) |>
    slice_head(n = 1) |>
    ungroup()

  unmapped <- setdiff(fg_genes, best_consolidated$SYMBOL)
  if (length(unmapped)) {
    best_consolidated <- bind_rows(
      best_consolidated,
      tibble(
        SYMBOL = unmapped, slim = "OTHER", consolidated = "Other",
        n_fg = NA_integer_, priority = 3L
      )
    )
  }

  small_cats <- best_consolidated |>
    count(consolidated) |>
    filter(n < min_cat_size, consolidated != "Other") |>
    pull(consolidated)
  if (length(small_cats)) {
    best_consolidated <- best_consolidated |>
      mutate(consolidated = ifelse(consolidated %in% small_cats, "Other", consolidated))
  }

  best_consolidated |>
    transmute(gene = SYMBOL, slim, consolidated) |>
    mutate(consolidated = factor(consolidated, levels = CONSOLIDATED_PATHWAY_ORDER))
}

# Cosine-blended ribbon polygon between two vertical bands (ported from YvO)
make_sigmoid_ribbon <- function(x0, x1, y0_top, y0_bot, y1_top, y1_bot,
                                n_pts = 50, ribbon_id) {
  t <- seq(0, 1, length.out = n_pts)
  blend <- (1 - cos(pi * t)) / 2
  tibble(
    x = c(x0 + (x1 - x0) * t, rev(x0 + (x1 - x0) * t)),
    y = c(
      y0_top + (y1_top - y0_top) * blend,
      rev(y0_bot + (y1_bot - y0_bot) * blend)
    ),
    ribbon_id = ribbon_id
  )
}

# Panel A: quadrant ORA scatter

# Per-protein HR-vs-LR logFC with pi-based significance class. Interaction-sig
# proteins are Divergent regardless of quadrant.
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
    distinct(gene, .keep_all = TRUE) |>
    mutate(
      quadrant = case_when(
        lfc_hi > 0 & lfc_lo > 0 ~ "Concordant Up",
        lfc_hi < 0 & lfc_lo < 0 ~ "Concordant Down",
        lfc_hi > 0 & lfc_lo < 0 ~ "HR Up / LR Down",
        TRUE ~ "HR Down / LR Up"
      ),
      sig_class = case_when(
        sig_int ~ "Divergent",
        sig_hi & sig_lo ~ "Sig Both",
        sig_hi ~ "Sig HR only",
        sig_lo ~ "Sig LR only",
        TRUE ~ "NS"
      ) |> factor(levels = names(SIG_COLORS_CONC))
    )
  q
}

QUAD_LEVELS_A <- c(
  "Concordant Up", "Concordant Down", "HR Up / LR Down", "HR Down / LR Up"
)

# Threshold-free per-quadrant ORA; padj_cutoff = 1 keeps leading terms even when
# a quadrant has no FDR-significant enrichment.
run_quadrant_ora <- function(quad_tbl, pw, n_show = 5) {
  universe <- quad_tbl$gene
  out <- lapply(QUAD_LEVELS_A, function(q) {
    genes <- quad_tbl$gene[quad_tbl$quadrant == q]
    if (length(genes) < 5) {
      return(NULL)
    }
    res <- tryCatch(
      run_ora_deduplicated(
        genes = genes, universe = universe, pathways = pw,
        jaccard_cutoff = 0.5, min_size = 15, max_size = 500, padj_cutoff = 1
      ),
      error = function(e) tibble()
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

# One quadrant's ORA as a half-bar block flanking its corner. Bars grow outward
# from the scatter edge; the pathway name sits just above (top) or below (bottom)
# each bar, anchored at the inner edge and read outward.
make_half_bars <- function(df, fill_color, side, ylim) {
  if (is.null(df) || nrow(df) == 0) {
    return(ggplot() +
      theme_void() +
      scale_y_continuous(limits = ylim, expand = c(0, 0)))
  }
  span <- max(abs(ylim))
  upper <- ylim[1] >= 0
  n_bars <- min(nrow(df), 5L)
  y_pos <- if (upper) {
    rev(seq(0.45, span - 0.4, length.out = 5))[seq_len(n_bars)]
  } else {
    seq(-0.45, -(span - 0.4), length.out = 5)[seq_len(n_bars)]
  }
  bars <- df |>
    arrange(desc(neg_log10_padj)) |>
    slice_head(n = 5) |>
    mutate(
      y = y_pos,
      bar_fill = ifelse(significant, scales::alpha(fill_color, 0.85),
        scales::alpha(fill_color, 0.3)
      ),
      star = sig_stars(padj)
    )
  x_max <- max(bars$neg_log10_padj)
  lab_hjust <- if (side == "left") 1 else 0
  lab_dy <- if (upper) 0.27 else -0.27
  star_nudge <- if (side == "left") -x_max * 0.04 else x_max * 0.04

  p <- ggplot(bars, aes(y = y)) +
    geom_rect(
      aes(xmin = 0, xmax = neg_log10_padj, ymin = y - 0.14, ymax = y + 0.14),
      fill = bars$bar_fill, color = "grey25", linewidth = 0.2
    ) +
    geom_text(aes(x = 0, y = y + lab_dy, label = pathway_label),
      hjust = lab_hjust, size = 2.5, fontface = "bold", color = "grey15"
    ) +
    geom_text(aes(x = neg_log10_padj + star_nudge, label = star),
      hjust = 0.5, size = 2.8, fontface = "bold", color = "black"
    ) +
    labs(x = expression(-log[10](p[adj])), y = NULL) +
    theme_minimal(base_size = 9) +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_blank(), axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = 7),
      axis.title.x = element_text(size = 7.5, face = "bold"),
      axis.line.x = element_line(color = "grey60", linewidth = 0.3)
    )
  x_scale <- if (side == "left") {
    scale_x_reverse(
      limits = c(x_max * 1.05, 0), breaks = scales::pretty_breaks(3),
      expand = expansion(mult = c(0.02, 0))
    )
  } else {
    scale_x_continuous(
      limits = c(0, x_max * 1.05), breaks = scales::pretty_breaks(3),
      expand = expansion(mult = c(0, 0.02))
    )
  }
  p + x_scale +
    scale_y_continuous(limits = ylim, expand = c(0, 0)) +
    coord_cartesian(clip = "off")
}

sig_key_A <- function() {
  kd <- tibble(
    label = c("Divergent", "Sig Both", "Sig HR only", "Sig LR only"),
    fill = unname(SIG_COLORS_CONC[c("Divergent", "Sig Both", "Sig HR only", "Sig LR only")]),
    x = c(1, 2, 3.1, 4.2)
  )
  ggplot(kd, aes(x, 0)) +
    geom_point(aes(fill = fill),
      shape = 21, size = 3, color = "grey45", stroke = 0.4
    ) +
    geom_text(aes(label = label), nudge_x = 0.08, hjust = 0, size = 2.7, color = "grey20") +
    scale_fill_identity() +
    scale_x_continuous(limits = c(0.5, 5.6), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-0.2, 0.2), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    theme_void()
}

# HR-vs-LR logFC scatter with flanking per-quadrant ORA half-bars; concordant
# diagonal tinted red, discordant off-diagonal blue; interaction-sig = Divergent.
panel_quadrant_ora <- function(quad_tbl, ora_df, cfg) {
  rho <- cor(quad_tbl$lfc_hi, quad_tbl$lfc_lo, method = "spearman", use = "complete.obs")
  sig_sub <- quad_tbl |> filter(sig_class != "NS")
  r_pear <- if (nrow(sig_sub) >= 3) {
    cor(sig_sub$lfc_hi, sig_sub$lfc_lo, use = "complete.obs")
  } else {
    NA_real_
  }
  lim <- max(abs(c(quad_tbl$lfc_hi, quad_tbl$lfc_lo))) * 1.05

  cnt <- quad_tbl |> count(quadrant, name = "total")
  sg <- sig_sub |> count(quadrant, name = "sig")
  qc <- left_join(cnt, sg, by = "quadrant") |>
    mutate(sig = tidyr::replace_na(sig, 0L))
  qlab <- function(q) {
    if (!q %in% qc$quadrant) {
      return("0/0")
    }
    sprintf("%d/%d", qc$sig[qc$quadrant == q], qc$total[qc$quadrant == q])
  }
  lab_df <- quad_tbl |>
    filter(sig_class != "NS") |>
    group_by(sig_class) |>
    slice_max(abs(lfc_hi) + abs(lfc_lo), n = 5, with_ties = FALSE) |>
    ungroup() |>
    mutate(label_fill = SIG_LABEL_FILL_CONC[as.character(sig_class)])

  scatter <- ggplot(quad_tbl, aes(lfc_hi, lfc_lo)) +
    annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf, fill = CONC_TINT, alpha = 0.6) +
    annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0, fill = CONC_TINT, alpha = 0.6) +
    annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0, fill = DISC_TINT, alpha = 0.6) +
    annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf, fill = DISC_TINT, alpha = 0.6) +
    geom_hline(yintercept = 0, color = "grey55", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "grey55", linewidth = 0.3) +
    geom_abline(slope = 1, linetype = "dashed", color = "black", linewidth = 0.3) +
    geom_point(data = ~ filter(.x, sig_class == "NS"), color = "grey80", size = 0.6, alpha = 0.35) +
    geom_point(
      data = ~ filter(.x, sig_class != "NS"), aes(fill = sig_class),
      shape = 21, size = 1.5, alpha = 0.9, color = "grey35", stroke = 0.25
    ) +
    scale_fill_manual(values = SIG_COLORS_CONC, guide = "none") +
    ggrepel::geom_label_repel(
      data = lab_df, aes(label = gene),
      fill = lab_df$label_fill, color = "white", size = 2.3, fontface = "italic",
      max.overlaps = 30, segment.size = 0.2, min.segment.length = 0,
      label.padding = unit(1, "pt"), label.r = unit(0.5, "pt"),
      box.padding = 0.3, seed = 42
    ) +
    annotate("label",
      x = lim, y = lim, label = sprintf("Concordant Up\n%s", qlab("Concordant Up")),
      hjust = 1, vjust = 1, size = 2.7, fontface = "bold", color = COMP_RED,
      fill = alpha("white", 0.9), label.padding = unit(2, "pt"), lineheight = 0.9
    ) +
    annotate("label",
      x = -lim, y = -lim, label = sprintf("%s\nConcordant Down", qlab("Concordant Down")),
      hjust = 0, vjust = 0, size = 2.7, fontface = "bold", color = COMP_RED,
      fill = alpha("white", 0.9), label.padding = unit(2, "pt"), lineheight = 0.9
    ) +
    annotate("label",
      x = -lim, y = lim, label = sprintf("HR Down / LR Up\n%s", qlab("HR Down / LR Up")),
      hjust = 0, vjust = 1, size = 2.7, fontface = "bold", color = COMP_BLUE,
      fill = alpha("white", 0.9), label.padding = unit(2, "pt"), lineheight = 0.9
    ) +
    annotate("label",
      x = lim, y = -lim, label = sprintf("%s\nHR Up / LR Down", qlab("HR Up / LR Down")),
      hjust = 1, vjust = 0, size = 2.7, fontface = "bold", color = COMP_BLUE,
      fill = alpha("white", 0.9), label.padding = unit(2, "pt"), lineheight = 0.9
    ) +
    coord_cartesian(xlim = c(-lim, lim), ylim = c(-lim, lim), expand = FALSE) +
    labs(x = cfg$labels$x, y = cfg$labels$y) +
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
  plot <- nw + scatter + ne + sw + se + sig_key_A() +
    plot_layout(design = design, widths = c(70, 105, 70), heights = c(85, 85, 12))
  list(
    plot = plot, rho = rho, r_pear = r_pear,
    n_total = nrow(quad_tbl), n_sig = nrow(sig_sub),
    n_enrich = sum(ora_df$significant, na.rm = TRUE)
  )
}

# Panel B: pattern heatmap + Sankey to consolidated GO-Slim + stacked count bars

panel_pattern_heatmap <- function(dep, cfg) {
  cx <- cfg$c_hi
  cy <- cfg$c_lo
  lfc_x_col <- paste0("logFC_", cx)
  lfc_y_col <- paste0("logFC_", cy)

  sig_hi_col <- paste0("sig_pi_", cx)
  sig_lo_col <- paste0("sig_pi_", cy)
  sig_int_col <- paste0("sig_pi_", cfg$c_int)
  sig_df <- dep |>
    transmute(
      gene,
      lfc_x = .data[[lfc_x_col]], lfc_y = .data[[lfc_y_col]],
      sig_hi = .data[[sig_hi_col]] != 0,
      sig_lo = .data[[sig_lo_col]] != 0,
      sig_int = .data[[sig_int_col]] != 0
    ) |>
    filter(!is.na(lfc_x), !is.na(lfc_y), sig_hi | sig_lo | sig_int) |>
    mutate(
      quadrant = case_when(
        lfc_x > 0 & lfc_y > 0 ~ "Concordant Up",
        lfc_x < 0 & lfc_y < 0 ~ "Concordant Down",
        TRUE ~ "Discordant"
      ),
      sig_cat = case_when(
        sig_hi & sig_lo ~ "Both",
        sig_hi ~ "HR",
        sig_lo ~ "LR",
        TRUE ~ "Inter."
      )
    ) |>
    distinct(gene, .keep_all = TRUE)

  go_result <- assign_go_slim_consolidated(sig_df$gene, dep$gene)
  sig_df <- sig_df |>
    left_join(go_result |> select(gene, consolidated), by = "gene") |>
    mutate(pathway = ifelse(is.na(consolidated), "Other", as.character(consolidated))) |>
    mutate(quadrant = factor(quadrant, levels = QUAD_ORDER_B)) |>
    arrange(quadrant, pathway, desc(lfc_x))

  n_total <- nrow(sig_df)

  ROW_H <- 0.078
  quad_counts <- sig_df |>
    count(quadrant, .drop = FALSE) |>
    mutate(quadrant = factor(quadrant, levels = QUAD_ORDER_B)) |>
    arrange(quadrant)

  y_pos <- numeric(n_total)
  quad_starts <- numeric(nrow(quad_counts))
  quad_ends <- numeric(nrow(quad_counts))
  idx <- 1
  current_y <- 0
  for (q in seq_len(nrow(quad_counts))) {
    nq <- quad_counts$n[q]
    quad_starts[q] <- current_y
    if (nq > 0) {
      for (p in seq_len(nq)) {
        y_pos[idx] <- current_y + (p - 0.5) * ROW_H
        idx <- idx + 1
      }
    }
    quad_ends[q] <- current_y + nq * ROW_H
    current_y <- current_y + nq * ROW_H
  }
  total_h <- current_y
  sig_df$y <- y_pos
  names(quad_starts) <- QUAD_ORDER_B
  names(quad_ends) <- QUAD_ORDER_B

  BAR_YMIN <- 0
  BAR_YMAX <- total_h

  pw_counts <- sig_df |>
    filter(pathway != "Other") |>
    count(pathway, name = "n_prot") |>
    arrange(desc(n_prot)) |>
    filter(n_prot >= 2)
  n_pw <- nrow(pw_counts)
  if (n_pw == 0) {
    return(list(
      plot = ggplot() +
        theme_void(), n_total = n_total, n_pw = 0L,
      data = sig_df, bar_data = tibble(), flow = tibble()
    ))
  }

  row_height <- (BAR_YMAX - BAR_YMIN) / n_pw
  pw_counts$y_center <- BAR_YMIN + row_height * (seq_len(n_pw) - 0.5)
  pw_counts$y_top <- BAR_YMIN + row_height * (seq_len(n_pw) - 1)
  pw_counts$y_bot <- BAR_YMIN + row_height * seq_len(n_pw)
  BAR_H <- row_height * 0.78

  dom_quad <- sig_df |>
    filter(pathway %in% pw_counts$pathway) |>
    group_by(pathway, quadrant) |>
    summarise(
      n = n(),
      lfc_sum = sum(abs(lfc_x) + abs(lfc_y)),
      .groups = "drop"
    ) |>
    group_by(pathway) |>
    arrange(pathway, desc(n), desc(lfc_sum)) |>
    slice_head(n = 1) |>
    ungroup() |>
    mutate(dom_quad = as.character(quadrant)) |>
    select(pathway, dom_quad)
  pw_counts <- pw_counts |> left_join(dom_quad, by = "pathway")

  STRIP_W <- 0.10
  TILE_W <- 0.80
  X_SIG <- 0.27
  X_COL1 <- X_SIG + STRIP_W / 2 + TILE_W / 2 + 0.01
  X_COL2 <- X_COL1 + TILE_W + 0.01
  X_QUAD <- X_COL2 + TILE_W / 2 + STRIP_W / 2 + 0.01
  HEAT_LEFT <- X_SIG - STRIP_W / 2
  HEAT_RIGHT <- X_QUAD + STRIP_W / 2

  X_SANK_L <- HEAT_RIGHT + 0.08
  X_SANK_R <- 3.2
  X_BAR_L <- 3.3
  BAR_SCALE <- 0.46
  count_max <- max(pw_counts$n_prot)
  X_BAR_MAX <- max(X_BAR_L + 14 * BAR_SCALE, X_BAR_L + count_max * BAR_SCALE)

  bar_data <- sig_df |>
    filter(pathway %in% pw_counts$pathway) |>
    count(pathway, quadrant, name = "n_seg") |>
    left_join(pw_counts |> select(pathway, y_center, n_prot), by = "pathway") |>
    group_by(pathway) |>
    arrange(pathway, desc(n_seg)) |>
    mutate(
      cum_n = cumsum(n_seg) - n_seg,
      xmin = X_BAR_L + cum_n * BAR_SCALE,
      xmax = X_BAR_L + (cum_n + n_seg) * BAR_SCALE,
      ymin = y_center - BAR_H / 2,
      ymax = y_center + BAR_H / 2
    ) |>
    ungroup()

  bg_stripes <- pw_counts |>
    transmute(
      xmin = X_BAR_L - 0.05, xmax = X_BAR_MAX + 2.0,
      ymin = y_top, ymax = y_bot,
      fill = QUAD_BG_B[dom_quad]
    )

  pw_labels <- pw_counts |>
    transmute(x = X_BAR_L + n_prot * BAR_SCALE + 0.08, y = y_center, label = pathway)

  count_ticks <- tibble(
    val = pretty(c(0, count_max), n = 4),
    x = X_BAR_L + val * BAR_SCALE,
    y_tick_top = BAR_YMAX, y_tick_bot = BAR_YMAX + ROW_H * 1.6,
    y_label = BAR_YMAX + ROW_H * 2.6
  ) |>
    filter(val >= 0, val <= count_max)

  flow_df <- sig_df |>
    filter(pathway %in% pw_counts$pathway) |>
    count(quadrant, pathway, name = "n_flow") |>
    filter(n_flow > 0)

  source_bands <- flow_df |>
    group_by(quadrant) |>
    mutate(
      total_q = sum(n_flow), frac = n_flow / total_q,
      q_start = quad_starts[as.character(quadrant)],
      q_end = quad_ends[as.character(quadrant)],
      q_height = q_end - q_start
    ) |>
    arrange(quadrant, match(pathway, pw_counts$pathway)) |>
    mutate(
      cum_frac = cumsum(frac) - frac,
      src_top = q_start + cum_frac * q_height,
      src_bot = q_start + (cum_frac + frac) * q_height
    ) |>
    ungroup()

  target_bands <- bar_data |>
    group_by(pathway) |>
    arrange(pathway, desc(n_seg)) |>
    mutate(
      frac = n_seg / sum(n_seg),
      cum_frac = cumsum(frac) - frac,
      tgt_top = ymin + cum_frac * (ymax - ymin),
      tgt_bot = ymin + (cum_frac + frac) * (ymax - ymin)
    ) |>
    ungroup() |>
    select(pathway, quadrant, tgt_top, tgt_bot)

  ribbon_df <- source_bands |>
    select(quadrant, pathway, n_flow, src_top, src_bot) |>
    left_join(target_bands, by = c("quadrant", "pathway"))

  all_ribbons <- purrr::pmap_dfr(ribbon_df, function(quadrant, pathway, n_flow,
                                                     src_top, src_bot, tgt_top, tgt_bot) {
    rid <- paste(quadrant, pathway, sep = "___")
    df <- make_sigmoid_ribbon(X_SANK_L, X_SANK_R, src_top, src_bot, tgt_top, tgt_bot,
      n_pts = 60, ribbon_id = rid
    )
    df$quadrant <- quadrant
    df$pathway <- pathway
    df
  })

  endpoint_bars <- bar_data |>
    transmute(
      xmin = X_SANK_R - 0.04, xmax = X_SANK_R + 0.04,
      ymin, ymax, quadrant = as.character(quadrant)
    )

  fc_max <- max(abs(c(sig_df$lfc_x, sig_df$lfc_y)), na.rm = TRUE)
  lfc_to_color <- function(v, fc_max) {
    v <- pmax(-fc_max, pmin(fc_max, v))
    ifelse(v >= 0,
      scales::seq_gradient_pal("#FFFFFF", unname(GROUP_COLORS["LR"]))(v / fc_max),
      scales::seq_gradient_pal(unname(GROUP_COLORS["HR"]), "#FFFFFF")((v + fc_max) / fc_max)
    )
  }

  heat_tiles <- bind_rows(
    sig_df |> transmute(x = X_COL1, y, w = TILE_W, h = ROW_H, fill = lfc_to_color(lfc_x, fc_max)),
    sig_df |> transmute(x = X_COL2, y, w = TILE_W, h = ROW_H, fill = lfc_to_color(lfc_y, fc_max))
  )
  sig_tiles <- sig_df |>
    transmute(x = X_SIG, y, w = STRIP_W, h = ROW_H, fill = SIG_COLORS_B[sig_cat])
  quad_tiles <- sig_df |>
    transmute(x = X_QUAD, y, w = STRIP_W, h = ROW_H, fill = QUAD_COLORS_B[as.character(quadrant)])

  divider_ys <- quad_ends[seq_len(length(QUAD_ORDER_B) - 1)]
  divider_ys <- divider_ys[divider_ys > 0 & divider_ys < total_h]

  col_headers <- tibble(
    x = c(X_COL1, X_COL2), y = total_h + ROW_H * 2.2,
    label = c(cfg$labels$hi, cfg$labels$lo),
    color = unname(c(CONTRAST_COLORS[cx], CONTRAST_COLORS[cy]))
  )

  n_g <- 50
  HEAT_MID <- (HEAT_LEFT + HEAT_RIGHT) / 2
  GRAD_HALFW <- (HEAT_RIGHT - HEAT_LEFT) * 0.30
  grad_xs <- seq(HEAT_MID - GRAD_HALFW, HEAT_MID + GRAD_HALFW, length.out = n_g)
  grad_h_legend <- tibble(
    xmin = grad_xs,
    xmax = lead(grad_xs, default = max(grad_xs) + diff(grad_xs)[1]),
    fill = lfc_to_color(seq(-fc_max, fc_max, length.out = n_g), fc_max)
  )
  GRAD_Y <- total_h + ROW_H * 3.0

  bar_mid <- (X_BAR_L + X_BAR_MAX) / 2
  KEY_Y_BASE <- BAR_YMAX + ROW_H * 15.5
  KEY_DY <- ROW_H * 3.8
  sig_cats <- c("HR", "LR", "Both", "Inter.")
  sig_cat_labels <- c(
    sprintf("Sig %s", cfg$labels$hi), sprintf("Sig %s", cfg$labels$lo),
    "Sig Both", "Interaction"
  )
  sig_key_df <- tibble(
    x = bar_mid - 1.43,
    y = KEY_Y_BASE + KEY_DY * seq_along(sig_cats),
    label = sig_cat_labels, fill = SIG_COLORS_B[sig_cats]
  )
  quad_key_df <- tibble(
    x = bar_mid + 0.6,
    y = KEY_Y_BASE + KEY_DY * seq_along(QUAD_ORDER_B),
    label = QUAD_ORDER_B, fill = QUAD_COLORS_B[QUAD_ORDER_B]
  )

  pB <- ggplot() +
    geom_rect(
      data = bg_stripes, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = bg_stripes$fill, color = "grey70", linewidth = 0.2
    ) +
    geom_rect(
      data = heat_tiles,
      aes(xmin = x - w / 2, xmax = x + w / 2, ymin = y - h / 2, ymax = y + h / 2),
      fill = heat_tiles$fill, color = NA
    ) +
    geom_rect(
      data = sig_tiles,
      aes(xmin = x - w / 2, xmax = x + w / 2, ymin = y - h / 2, ymax = y + h / 2),
      fill = sig_tiles$fill, color = NA
    ) +
    geom_rect(
      data = quad_tiles,
      aes(xmin = x - w / 2, xmax = x + w / 2, ymin = y - h / 2, ymax = y + h / 2),
      fill = quad_tiles$fill, color = NA
    ) +
    geom_segment(
      data = tibble(y = divider_ys),
      aes(x = X_SIG - STRIP_W / 2, xend = X_QUAD + STRIP_W / 2, y = y, yend = y),
      color = "grey30", linewidth = 0.4
    ) +
    geom_text(
      data = col_headers, aes(x = x, y = y, label = label),
      size = 3.4, fontface = "bold", color = col_headers$color
    ) +
    geom_polygon(
      data = all_ribbons, aes(x = x, y = y, group = ribbon_id),
      fill = QUAD_COLORS_B[all_ribbons$quadrant], alpha = 0.40, color = NA
    ) +
    geom_rect(
      data = endpoint_bars, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = ENDPOINT_COLORS_B[endpoint_bars$quadrant], color = NA
    ) +
    geom_rect(
      data = bar_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = QUAD_COLORS_B[as.character(bar_data$quadrant)], color = "black", linewidth = 0.3
    ) +
    geom_text(
      data = bar_data, aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = n_seg),
      size = 2.6, fontface = "bold", color = "white"
    ) +
    geom_text(
      data = pw_labels, aes(x = x, y = y, label = label),
      size = 3.0, hjust = 0, fontface = "bold", color = "grey15", lineheight = 0.8
    ) +
    annotate("segment",
      x = X_BAR_L - 0.05, xend = X_BAR_MAX + 2.0,
      y = BAR_YMAX, yend = BAR_YMAX, color = "grey20", linewidth = 0.5
    ) +
    geom_segment(
      data = count_ticks, aes(x = x, xend = x, y = y_tick_top, yend = y_tick_bot),
      color = "grey20", linewidth = 0.3
    ) +
    geom_text(
      data = count_ticks, aes(x = x, y = y_label, label = val),
      size = 2.7, fontface = "bold", color = "grey20"
    ) +
    annotate("text",
      x = X_BAR_L + 7.5 * BAR_SCALE, y = BAR_YMAX + ROW_H * 4.4,
      label = "Protein count", size = 3.2, fontface = "bold", color = "grey20"
    ) +
    geom_rect(
      data = grad_h_legend, aes(xmin = xmin, xmax = xmax, ymin = GRAD_Y, ymax = GRAD_Y + ROW_H * 1.4),
      fill = grad_h_legend$fill, color = NA
    ) +
    annotate("text",
      x = HEAT_MID, y = GRAD_Y + ROW_H * 3.0, label = "logFC", size = 3.0,
      fontface = "bold", color = "grey20"
    ) +
    annotate("text", x = HEAT_LEFT, y = GRAD_Y + ROW_H * 3.0, label = sprintf("%.1f", -fc_max), size = 2.6, color = "grey30") +
    annotate("text", x = HEAT_RIGHT, y = GRAD_Y + ROW_H * 3.0, label = sprintf("+%.1f", fc_max), size = 2.6, color = "grey30") +
    geom_point(
      data = sig_key_df, aes(x = x, y = y), shape = 22, size = 3.2,
      fill = sig_key_df$fill, color = "grey30", stroke = 0.3
    ) +
    geom_text(
      data = sig_key_df, aes(x = x + 0.12, y = y, label = label),
      hjust = 0, size = 2.7, color = "grey20"
    ) +
    geom_point(
      data = quad_key_df, aes(x = x, y = y), shape = 22, size = 3.2,
      fill = quad_key_df$fill, color = "grey30", stroke = 0.3
    ) +
    geom_text(
      data = quad_key_df, aes(x = x + 0.12, y = y, label = label),
      hjust = 0, size = 2.7, color = "grey20"
    ) +
    scale_y_reverse() +
    coord_cartesian(
      xlim = c(0.0, X_BAR_MAX + 2.0),
      ylim = c(BAR_YMAX + ROW_H * 20, -ROW_H * 0.05),
      expand = FALSE
    ) +
    theme_void() +
    theme(plot.margin = margin(2, 2, 2, 2, "mm"))

  list(plot = pB, n_total = n_total, n_pw = n_pw, data = sig_df, bar_data = bar_data, flow = flow_df)
}

# Panel C: fry rotation test

running_es <- function(t_vals, in_set) {
  n <- length(t_vals)
  n_h <- sum(in_set)
  if (n_h == 0) {
    return(rep(0, n))
  }
  hit_w <- ifelse(in_set, abs(t_vals), 0)
  cumsum(ifelse(in_set, hit_w / sum(hit_w), -1 / (n - n_h)))
}

# fry rotation test of the HR pi-sig up/down DEP sets along the LR moderated-t
# ranking. fry carries the within-subject block + duplicateCorrelation, so it is
# valid for this repeated-measures design. Sets are HR pi-sig; the test axis is
# LR_T2-LR_T1 (training) or LR_T3-LR_T2 (acute); block = Subject_ID.
run_fry_concordance <- function(da, dep, pw, c_hi, c_lo, lo_levels) {
  set.seed(42)
  mat <- as.matrix(da$data)
  meta <- da$metadata
  meta$Group_Time <- factor(meta$Group_Time)
  design <- model.matrix(~ 0 + Group_Time, data = meta)
  colnames(design) <- sub("^Group_Time", "", colnames(design))
  block <- meta$Subject_ID
  corfit <- duplicateCorrelation(mat, design, block = block)
  cor_within <- corfit$consensus.correlation
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
    block = block, correlation = cor_within
  )
  res <- tibble::rownames_to_column(as.data.frame(res), "set")
  res$expected <- res$set
  res$consistent <- res$Direction == res$expected
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
    r <- tryCatch(
      run_ora_deduplicated(genes, rk$gene, pw,
        jaccard_cutoff = 0.5, min_size = 15, max_size = 500, padj_cutoff = 1
      ),
      error = function(e) tibble()
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
    n_all = nrow(rk), cor_within = cor_within
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
    geom_text(aes(label = sig_stars(padj)), hjust = -0.2, size = 2.6) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(x = expression(-log[10](p[adj])), y = NULL) +
    FIG_THEME +
    theme(axis.text.y = element_text(size = 7))
}

# Running-ES curves + barcode rugs + ranking metric + leading-edge ORA in the YvO
# fry layout (ES + barcode per direction stacked, ORA flanking).
panel_fry <- function(fry_out, cfg) {
  rk <- fry_out$rank
  res <- fry_out$results
  n_all <- fry_out$n_all
  fr <- function(s) res[res$set == s, , drop = FALSE]

  es_plot <- function(es_col, color, title, set) {
    row <- fr(set)
    is_sig <- nrow(row) && !is.na(row$PValue) && row$PValue < 0.05
    col <- if (is_sig) color else scales::alpha(color, 0.4)
    lab <- if (nrow(row)) {
      sprintf(
        "fry %s, %s (n = %d)%s", row$Direction, fmt_p(row$PValue), row$n,
        if (row$consistent) "" else " \u2717"
      )
    } else {
      set
    }
    ggplot(rk, aes(rank, .data[[es_col]])) +
      geom_area(fill = scales::alpha(col, 0.15), color = NA) +
      geom_line(color = col, linewidth = 0.6) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.3) +
      annotate("text",
        x = n_all, y = Inf, label = lab, hjust = 1.05, vjust = 1.6,
        size = 2.6, fontface = "bold", color = "grey20"
      ) +
      labs(title = title, x = NULL, y = "ES") +
      scale_x_continuous(limits = c(1, n_all), expand = c(0.005, 0)) +
      FIG_THEME +
      theme(
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.title = element_text(size = 9.5, face = "bold"),
        plot.margin = margin(0, 1, 0, 0, "mm")
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
    labs(x = sprintf("Rank (%s moderated t, n = %d)", cfg$labels$lo, n_all), y = NULL) +
    scale_x_continuous(limits = c(1, n_all), expand = c(0.005, 0)) +
    FIG_THEME +
    theme(plot.margin = margin(0, 1, 1, 0, "mm"))

  design <- c(
    patchwork::area(1, 1), patchwork::area(2, 1), patchwork::area(3, 1),
    patchwork::area(4, 1), patchwork::area(5, 1),
    patchwork::area(1, 2, 2, 2), patchwork::area(3, 2, 4, 2)
  )
  es_plot("es_up", COMP_RED, sprintf("HR-up DEPs \u2192 %s rank", cfg$labels$lo), "Up") +
    bc_plot("in_up", COMP_RED, "Up") +
    es_plot("es_down", COMP_BLUE, sprintf("HR-down DEPs \u2192 %s rank", cfg$labels$lo), "Down") +
    bc_plot("in_down", COMP_BLUE, "Down") +
    rank_curve +
    ora_bars(fry_out$ora_up, COMP_RED) +
    ora_bars(fry_out$ora_down, COMP_BLUE) +
    plot_layout(
      design = design, heights = c(2, 0.25, 2, 0.25, 0.7), widths = c(1.7, 1.5)
    )
}

# Panel D: pathway NES scatter (Hallmark + GO Slim)

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
    ) |>
    mutate(set_size = coalesce(.data[[paste0("size_", c_hi)]], .data[[paste0("size_", c_lo)]]))
  int <- cache |>
    filter(contrast == c_int, database %in% keep_db) |>
    select(pathway, NES_int = NES, padj_int = padj)
  left_join(wide, int, by = "pathway")
}

# NES(HR) vs NES(LR); ref_slope +1, diagonal concordant. Shape = database, size =
# set size; interaction-divergent pathways (cached interaction padj < 0.05)
# highlighted. Spearman rho + Fisher-z CI in the subtitle.
panel_nes_scatter <- function(nes_wide, c_hi, c_lo, cfg) {
  nes_x <- paste0("NES_", c_hi)
  nes_y <- paste0("NES_", c_lo)
  x <- nes_wide[[nes_x]]
  y <- nes_wide[[nes_y]]
  ct <- suppressWarnings(cor.test(x, y, method = "spearman"))
  rho <- unname(ct$estimate)
  ci <- fisher_z_ci_spearman(rho, nrow(nes_wide))
  lim <- max(abs(c(x, y))) * 1.35

  d <- nes_wide |>
    mutate(
      nx = x, ny = y,
      divergent = !is.na(padj_int) & padj_int < 0.05,
      class = if_else(divergent, "Divergent", "Concordant")
    )
  n_conc <- sum(sign(d$nx) == sign(d$ny))
  conc_frac <- n_conc / nrow(d)

  lab_d <- d |>
    filter(divergent | (sign(nx) == sign(ny) & abs(nx) + abs(ny) > 2)) |>
    slice_max(abs(nx) + abs(ny), n = 14)

  tints <- tibble(
    xmin = c(0, -lim, 0, -lim), xmax = c(lim, 0, lim, 0),
    ymin = c(0, -lim, -lim, 0), ymax = c(lim, 0, 0, lim),
    fill = c(CONC_TINT, CONC_TINT, DISC_TINT, DISC_TINT)
  )

  p <- ggplot(d, aes(nx, ny)) +
    geom_rect(
      data = tints, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
      inherit.aes = FALSE, alpha = 0.5
    ) +
    scale_fill_identity() +
    geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "grey60", linewidth = 0.3) +
    geom_abline(slope = 1, linetype = "dashed", color = "black", linewidth = 0.3) +
    geom_point(
      data = ~ filter(.x, !divergent),
      aes(shape = database, size = set_size), fill = "grey65", color = "grey45",
      alpha = 0.7, stroke = 0.2
    ) +
    geom_point(
      data = ~ filter(.x, divergent),
      aes(shape = database, size = set_size), fill = "#E08214", color = "black",
      alpha = 0.9, stroke = 0.4
    ) +
    scale_shape_manual(values = c(Hallmark = 24, "GO Slim" = 21), name = "Database") +
    scale_size_continuous(range = c(1.5, 5), name = "Set size", breaks = c(20, 50, 100, 200)) +
    ggrepel::geom_label_repel(
      data = lab_d, aes(label = clean_pathway_name(pathway, 26)),
      color = "black", fill = alpha("white", 0.85), size = 2.2, fontface = "bold",
      lineheight = 0.82, max.overlaps = Inf, segment.size = 0.2, segment.color = "grey50",
      min.segment.length = 0, box.padding = 0.4, label.padding = unit(1, "pt"),
      label.r = unit(0.5, "pt"), label.size = 0.2, seed = 42
    ) +
    annotate("label",
      x = lim, y = lim, label = sprintf("Concordant Up  n = %d", sum(d$nx > 0 & d$ny > 0)),
      hjust = 1, vjust = 1, size = 2.6, fontface = "bold", color = COMP_RED,
      fill = alpha("white", 0.92), label.padding = unit(2, "pt")
    ) +
    annotate("label",
      x = -lim, y = -lim, label = sprintf("Concordant Down  n = %d", sum(d$nx < 0 & d$ny < 0)),
      hjust = 0, vjust = 0, size = 2.6, fontface = "bold", color = COMP_RED,
      fill = alpha("white", 0.92), label.padding = unit(2, "pt")
    ) +
    annotate("label",
      x = -lim, y = lim, label = sprintf("Discordant  n = %d", sum(d$nx < 0 & d$ny > 0)),
      hjust = 0, vjust = 1, size = 2.6, fontface = "bold", color = COMP_BLUE,
      fill = alpha("white", 0.92), label.padding = unit(2, "pt")
    ) +
    annotate("label",
      x = lim, y = -lim, label = sprintf("Discordant  n = %d", sum(d$nx > 0 & d$ny < 0)),
      hjust = 1, vjust = 0, size = 2.6, fontface = "bold", color = COMP_BLUE,
      fill = alpha("white", 0.92), label.padding = unit(2, "pt")
    ) +
    coord_fixed(ratio = 1, xlim = c(-lim, lim), ylim = c(-lim, lim)) +
    labs(
      x = sprintf("NES (%s)", cfg$labels$x_short),
      y = sprintf("NES (%s)", cfg$labels$y_short),
      subtitle = sprintf(
        "GO Slim + Hallmark – %d pathways – rho = %.2f [%.2f, %.2f] – %.0f%% concordant",
        nrow(d), rho, ci["lo"], ci["hi"], conc_frac * 100
      )
    ) +
    FIG_THEME +
    theme(
      legend.position = "bottom", legend.box = "horizontal",
      legend.key.size = unit(3, "mm")
    ) +
    guides(
      shape = guide_legend(order = 1, override.aes = list(size = 3, fill = "grey50")),
      size = guide_legend(order = 2)
    )
  list(plot = p, rho = rho, ci = ci, data = d, conc_frac = conc_frac)
}

# Panel E: RRHO2 threshold-free concordance

JET_COLORS <- c(
  "#00007F", "blue", "#007FFF", "cyan", "#7FFF7F",
  "yellow", "#FF7F00", "red", "#7F0000"
)

# RRHO2 of the two moderated-t rankings (t_HR vs t_LR). Warm corners = concordant
# gene regulation; UU (top-left) + DD (bottom-right) are the concordant hotspots.
panel_rrho2 <- function(dep, c_hi, c_lo, cfg) {
  rr <- tibble(
    gene = dep$gene,
    t1 = dep[[paste0("t_", c_hi)]],
    t2 = dep[[paste0("t_", c_lo)]]
  ) |>
    filter(!is.na(gene), !is.na(t1), !is.na(t2)) |>
    distinct(gene, .keep_all = TRUE)
  n_shared <- nrow(rr)

  obj <- RRHO2_initialize(
    data.frame(gene = rr$gene, score = rr$t1, stringsAsFactors = FALSE),
    data.frame(gene = rr$gene, score = rr$t2, stringsAsFactors = FALSE),
    labels = c(cfg$labels$hi, cfg$labels$lo), log10.ind = TRUE,
    multipleTesting = "none", boundary = 0.02, method = "hyper", stepsize = 20
  )
  hmat <- obj$hypermat
  nr <- nrow(hmat)
  nc <- ncol(hmat)

  hotspot <- list(
    UU = obj$genelist_uu$gene_list_overlap_uu,
    DD = obj$genelist_dd$gene_list_overlap_dd,
    UD = obj$genelist_ud$gene_list_overlap_ud,
    DU = obj$genelist_du$gene_list_overlap_du
  )
  n_UU <- length(hotspot$UU)
  n_DD <- length(hotspot$DD)

  hm_df <- expand.grid(row = seq_len(nr), col = seq_len(nc)) |>
    mutate(neg_log10_p = as.vector(hmat))
  max_val <- max(hm_df$neg_log10_p, na.rm = TRUE)

  ql <- list(
    UU = "Concordant Up", DD = "Concordant Down",
    UD = sprintf("%s\u2191 %s\u2193", cfg$labels$hi, cfg$labels$lo),
    DU = sprintf("%s\u2193 %s\u2191", cfg$labels$hi, cfg$labels$lo)
  )
  LABEL_FILL <- scales::alpha("white", 0.85)

  plot <- ggplot(hm_df, aes(x = row, y = col, fill = neg_log10_p)) +
    geom_raster() +
    scale_fill_gradientn(
      colors = JET_COLORS, limits = c(0, max_val), na.value = "white",
      name = expression(-log[10](P)),
      guide = guide_colorbar(
        barwidth = unit(20, "mm"), barheight = unit(2.5, "mm"),
        title.position = "left", title.vjust = 1
      )
    ) +
    annotate("label",
      x = 1, y = 1, label = ql$UU, color = "grey15", fill = LABEL_FILL,
      label.size = 0, fontface = "bold", size = 2.6, hjust = 0, vjust = 0
    ) +
    annotate("label",
      x = nr, y = nc, label = ql$DD, color = "grey15", fill = LABEL_FILL,
      label.size = 0, fontface = "bold", size = 2.6, hjust = 1, vjust = 1
    ) +
    annotate("label",
      x = 1, y = nc, label = ql$UD, color = "grey15", fill = LABEL_FILL,
      label.size = 0, fontface = "bold", size = 2.6, hjust = 0, vjust = 1
    ) +
    annotate("label",
      x = nr, y = 1, label = ql$DU, color = "grey15", fill = LABEL_FILL,
      label.size = 0, fontface = "bold", size = 2.6, hjust = 1, vjust = 0
    ) +
    scale_x_continuous(expand = expansion(mult = 0.015)) +
    scale_y_continuous(expand = expansion(mult = 0.015)) +
    labs(
      x = sprintf("%s rank (Up \u2192 Down)", cfg$labels$hi),
      y = sprintf("%s rank (Up \u2192 Down)", cfg$labels$lo)
    ) +
    FIG_THEME +
    theme(
      axis.text = element_blank(), axis.ticks = element_blank(),
      panel.border = element_blank(), panel.grid.major = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(size = 7)
    ) +
    coord_fixed(ratio = 1, clip = "off")

  hotspot_export <- bind_rows(
    tibble(quadrant = "UU", gene = hotspot$UU),
    tibble(quadrant = "DD", gene = hotspot$DD),
    tibble(quadrant = "UD", gene = hotspot$UD),
    tibble(quadrant = "DU", gene = hotspot$DU)
  )
  list(
    plot = plot, n_shared = n_shared, n_UU = n_UU, n_DD = n_DD,
    n_concordant = n_UU + n_DD, hotspot = hotspot_export
  )
}

# Supplements

# Quadrant split across significance thresholds.
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
    geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 2.6, color = "white") +
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

# Bootstrap CI on the HR-vs-LR logFC Spearman concordance (proteins resampled).
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
    annotate("rect", xmin = ci[1], xmax = ci[2], ymin = 0, ymax = Inf, alpha = 0.15, fill = unname(GROUP_COLORS["HR"])) +
    geom_vline(xintercept = obs, linetype = "dashed", color = unname(GROUP_COLORS["LR"]), linewidth = 0.6) +
    labs(
      x = "Spearman rho (logFC HR vs LR)", y = "Bootstrap reps",
      subtitle = sprintf("rho = %.2f [%.2f, %.2f] (95%% CI, %d protein resamples)", obs, ci[1], ci[2], reps)
    ) +
    FIG_THEME
  list(plot = plot, data = tibble(
    stat = "spearman_rho", observed = obs, ci_lo = ci[1], ci_hi = ci[2], reps = reps
  ))
}

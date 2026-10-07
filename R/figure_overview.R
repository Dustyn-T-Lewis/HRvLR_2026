# Figure 1 panels for any arc that compares two groups (HR/LR, Hi/Lo, C1/C2). Each builder takes
# saved outputs and the two group names, and returns list(plot, data). figure_theme.R and
# figure_panels.R must be sourced first, with dplyr, tidyr, purrr, tibble and ggplot2 attached.

# The nine contrasts of a two-group design with their panel family and row label: the group for a
# within-group window, the time point for a group difference, the window for an interaction.
overview_contrasts <- function(groups) {
  between_family <- paste(groups, collapse = " - ")
  versus <- paste(groups, collapse = "v")
  tibble(
    contrast = c(
      paste0("Training_", groups), paste0("Acute_", groups),
      paste0(c("Baseline_", "Trained_", "Acute_"), versus),
      "Training_Interaction", "Acute_Interaction"
    ),
    group = c(
      "Training", "Training", "Acute", "Acute", rep(between_family, 3),
      "Interaction", "Interaction"
    ),
    label = c(groups, groups, "T1", "T2", "T3", "Training", "Acute")
  ) |>
    mutate(
      group = factor(group, c("Training", "Acute", between_family, "Interaction")),
      contrast = factor(contrast, rev(contrast))
    )
}

# Proteins per contrast. The light bar counts every protein at p < 0.05 (or Pi < 0.05); the dark
# bar, drawn over it from zero, counts the Pi < 0.05 subset. Each row's name sits in a box in its
# family's tint and border. validity, named by contrast, adds a tag to each row label.
panel_protein_counts <- function(proteins, groups, n_tested, validity = NULL) {
  contrast_labels <- overview_contrasts(groups)
  if (!is.null(validity)) {
    tag <- c(circular = " (circ.)", caution = " (caut.)", clean = "")
    contrast_labels <- mutate(contrast_labels,
      label = paste0(label, tag[validity[as.character(contrast)]])
    )
  }
  counts <- proteins |>
    mutate(direction = if_else(logFC > 0, "Up", "Down")) |>
    summarise(
      pi = sum(pi_score < 0.05), nominal = sum(p < 0.05 | pi_score < 0.05),
      fdr = sum(fdr < 0.05), .by = c(contrast, direction)
    ) |>
    complete(contrast, direction, fill = list(pi = 0, nominal = 0, fdr = 0)) |>
    inner_join(mutate(contrast_labels, contrast = as.character(contrast)), by = "contrast")
  pi_text <- "\u03a0 < 0.05"
  bars <- counts |>
    pivot_longer(c(nominal, pi), names_to = "tier", values_to = "n") |>
    mutate(
      signed = if_else(direction == "Up", n, -n),
      contrast = factor(contrast, levels(contrast_labels$contrast)),
      fill = paste(direction, if_else(tier == "pi", pi_text, "p < 0.05"))
    ) |>
    # The Pi bar is drawn last, over the nominal bar.
    arrange(tier == "pi")
  bar_limit <- max(counts$nominal)
  bar_fills <- set_names(
    c(direction_colours[["Down"]], "#F1B6DA", "#B8E186", direction_colours[["Up"]]),
    paste(c("Down", "Down", "Up", "Up"), c(pi_text, "p < 0.05", "p < 0.05", pi_text))
  )
  row_boxes <- mutate(contrast_labels, family = as.character(group), x = -bar_limit * 1.06)
  n_fdr <- sum(counts$fdr)
  plot <- ggplot(bars, aes(signed, contrast)) +
    geom_col(aes(fill = fill), width = 0.7, position = "identity") +
    geom_vline(xintercept = 0, linewidth = 0.25) +
    # Families share the fill scale with the bars; the legend lists the bars only.
    geom_label(
      aes(x, contrast, label = label, fill = family, colour = family), row_boxes,
      hjust = 1, size = 2.1, fontface = "bold", text.colour = "black",
      linewidth = 0.35, label.r = unit(0.6, "mm"), label.padding = unit(0.7, "mm")
    ) +
    scale_fill_manual(
      values = c(bar_fills, contrast_tone("fill")), breaks = names(bar_fills), name = NULL
    ) +
    scale_colour_manual(values = contrast_tone("line"), guide = "none") +
    scale_x_continuous(
      labels = abs, limits = c(-bar_limit * 1.45, bar_limit), breaks = seq(-100, 100, 50)
    ) +
    facet_grid(group ~ ., scales = "free_y", space = "free_y") +
    coord_cartesian(clip = "off") +
    labs(
      x = "proteins (down | up)", y = NULL, title = "Proteins per contrast",
      subtitle = sprintf(
        "%s tested; %s at FDR < 0.05", scales::comma(n_tested),
        if (n_fdr == 0) "none" else scales::comma(n_fdr)
      )
    ) +
    figure_theme +
    theme(
      strip.text = element_blank(), axis.text.y = element_blank(),
      legend.position = "bottom", panel.grid.major.y = element_blank(),
      panel.spacing.y = unit(1.2, "mm")
    ) +
    guides(fill = guide_legend(nrow = 2))
  list(plot = plot, data = select(counts, group, label, contrast, direction, pi, nominal, fdr))
}

# Proteome and phenotype PCA side by side, each with its PERMANOVA (from the protein-contrasts
# stage). The phenotype is the six size and strength traits at T1 and T2, scaled; a biopsy missing
# a trait is left out rather than filled. Groups are not tested on the phenotype when they were
# split on it. group_of relabels people for a cluster arc.
panel_pca_pair <- function(pca, variance, permanova, phenotype, groups, group_of = NULL,
                           n_proteins = NULL) {
  shapes <- set_names(unname(arm_shapes), groups)
  test_row <- \(data_, term_, window_) {
    filter(permanova, data == data_, term == term_, window == window_, arms == "both")
  }
  permanova_text <- \(label, test) sprintf("%s R\u00b2 %.2f, p = %.2g", label, test$r2, test$p)
  time_test <- test_row("proteome", "time", "all")
  arm_test <- test_row("proteome", "arm", "subject mean")
  trait_time_test <- test_row("phenotype", "time", "all")
  proteome <- pca_panel(
    pca, variance,
    "Proteome", sprintf("%d biopsies, %s proteins", nrow(pca), scales::comma(n_proteins)),
    paste(
      "PERMANOVA", permanova_text("time", time_test), permanova_text("group", arm_test),
      sep = "\n"
    ),
    shapes = shapes
  )
  tables <- phenotype_tables(phenotype)
  traits <- filter(tables$biopsy, if_all(all_of(tables$traits), \(x) !is.na(x)))
  if (!is.null(group_of)) {
    traits <- mutate(traits, arm = unname(group_of[subject])) |> filter(!is.na(arm))
  }
  trait_pca <- prcomp(scale(traits[tables$traits]))
  scores <- bind_cols(select(traits, subject, arm, timepoint), as_tibble(trait_pca$x[, 1:2]))
  phenotype_panel <- pca_panel(
    scores, summary(trait_pca)$importance[2, ],
    "Phenotype", sprintf("%d biopsies, six size and strength traits", nrow(traits)),
    paste("PERMANOVA", permanova_text("time", trait_time_test), sep = "\n"),
    shapes = shapes
  ) +
    # No T3 biopsy has traits, so its key would differ from the proteome's and print twice.
    guides(colour = "none", fill = "none", shape = "none")
  plot <- (proteome | phenotype_panel) +
    patchwork::plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  list(
    plot = plot, proteome = pca, phenotype = scores,
    tests = list(time = time_test, arm = arm_test, trait_time = trait_time_test)
  )
}

# Pathways at fgsea FDR < 0.05 in both groups of a window (group difference: at least two of the
# three time points; interaction: either interaction, strongest five), from four collections, each
# shown in all nine contrasts. A set sharing half its measured proteins with a stronger set already
# shown is dropped. Opposite signs count only among the significant columns.
panel_pathway_grid <- function(set_tests, gene_sets, groups) {
  between_family <- paste(groups, collapse = " - ")
  versus <- paste(groups, collapse = "v")
  windows <- list(
    Training = set_names(paste0("Training_", groups), groups),
    Acute = set_names(paste0("Acute_", groups), groups),
    set_names(paste0(c("Baseline_", "Trained_", "Acute_"), versus), c("T1", "T2", "T3")),
    Interaction = c(Training = "Training_Interaction", Acute = "Acute_Interaction")
  )
  names(windows)[3] <- between_family
  needed <- set_names(c(2, 2, 2, 1), names(windows))
  row_cap <- set_names(c(Inf, Inf, Inf, 5), names(windows))
  collections <- c(
    Hallmark = "#332288", Reactome = "#88CCEE", KEGG_Legacy = "#DDCC77", GOBP = "#777777"
  )
  members <- gene_sets$sets
  fgsea <- set_tests$set_tests |>
    filter(method == "fgsea", collection %in% names(collections)) |>
    select(contrast, set_id, collection, pathway, nes, fdr, leading_edge)
  select_window <- function(name) {
    rows <- filter(fgsea, contrast %in% windows[[name]])
    eligible <- rows |>
      summarise(
        hits = sum(fdr < 0.05),
        same_sign = n_distinct(sign(nes[fdr < 0.05])) == 1,
        strength = sum(abs(nes)),
        # Several Reactome signalling sets carry only proteasome subunits; the label says so.
        proteasome = mean(startsWith(
          as.character(unique(unlist(leading_edge[fdr < 0.05]))), "PSM"
        )),
        .by = c(set_id, pathway)
      ) |>
      filter(hits >= needed[[name]]) |>
      arrange(desc(same_sign), desc(strength))
    kept <- head(prune_overlapping(eligible$set_id, members), row_cap[[name]])
    eligible |>
      filter(set_id %in% kept) |>
      transmute(
        facet = name, set_id, order = match(set_id, kept),
        label = paste0(
          if_else(proteasome >= 0.5, "Proteasome: ", ""), pathway_label(pathway),
          if_else(same_sign, "", " *")
        )
      )
  }
  columns <- imap(windows, \(contrasts, window) {
    tibble(window, arm = names(contrasts), contrast = unname(contrasts))
  }) |>
    list_rbind() |>
    mutate(
      column = factor(paste(window, arm), c(" ", paste(window, arm))),
      markup = in_colour(arm, arm_shade(window, arm, groups)),
      window = factor(window, c("", names(windows)))
    )
  pathways <- map(names(windows), select_window) |>
    list_rbind() |>
    inner_join(fgsea, by = "set_id", relationship = "many-to-many") |>
    inner_join(columns, by = "contrast") |>
    mutate(
      facet = factor(facet, names(windows)), row = paste(facet, set_id),
      log_fdr = pmin(-log10(fdr), 12)
    )
  # A narrow first column marks each set's database.
  databases <- distinct(pathways, row, facet, order, collection) |>
    mutate(
      window = factor("", levels(columns$window)), column = factor(" ", levels(columns$column))
    )
  nes_limit <- ceiling(max(abs(pathways$nes)) * 10) / 10
  plot <- pathways |>
    mutate(row = reorder(row, -order)) |>
    ggplot(aes(column, row)) +
    geom_point(aes(size = log_fdr, fill = nes), shape = 21, stroke = 0.2) +
    geom_point(
      aes(colour = collection), mutate(databases, row = reorder(row, -order)),
      shape = 15, size = 1.8
    ) +
    scale_colour_manual(
      values = collections, breaks = names(collections),
      labels = c("Hallmark", "Reactome", "KEGG", "GO BP"), name = "Database"
    ) +
    scale_y_discrete(labels = deframe(distinct(pathways, row, label))) +
    scale_x_discrete(labels = c(` ` = "", deframe(select(columns, column, markup)))) +
    scale_fill_distiller(
      palette = "PiYG", direction = 1, limits = c(-nes_limit, nes_limit), name = "NES"
    ) +
    scale_size_area(
      max_size = 2.8, limits = c(0, 12), breaks = c(2, 6, 12),
      labels = c("2", "6", "12+"), name = "-log10 FDR"
    ) +
    guides(
      fill = guide_colourbar(order = 1), size = guide_legend(order = 2),
      colour = guide_legend(order = 3)
    ) +
    ggh4x::facet_grid2(
      facet ~ window,
      scales = "free", space = "free", switch = "y",
      strip = contrast_strips(x = levels(columns$window), y = names(windows))
    ) +
    labs(
      x = NULL, y = NULL, title = "Pathways shared within each window",
      subtitle = paste(
        "fgsea FDR < 0.05 in both groups of a window (interaction: either, top 5),",
        "shown in all nine contrasts; * opposite signs"
      )
    ) +
    figure_theme +
    theme(
      axis.text.y = element_text(size = 5.2),
      axis.text.x.bottom = ggtext::element_markdown(),
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 90, size = 5.5, face = "bold"),
      strip.text.x = element_text(size = 6, face = "bold"),
      panel.spacing = unit(1, "mm")
    )
  list(
    plot = plot,
    data = transmute(
      pathways,
      selected_in = as.character(facet), window = as.character(window),
      column = as.character(arm), set_id, collection, label, nes, fdr
    ),
    n_columns = nrow(columns)
  )
}

# Figure 1 for a two-group arc: the arc's own panel a (design or clustering), then b proteins per
# contrast, c proteome and phenotype PCA, d the pathway grid, and e-f each group's training and
# acute sets ranked in the other group. inputs must name contrasts (the protein-contrasts
# workbook), set_tests, gene_sets and phenotype. Every plotted number is checked against its
# source table before anything is written.
figure_overview <- function(panel_a, inputs, paths, groups, out_dir, stem, a_sheets = list(),
                            a_holds = character(), note = NULL, validity = NULL,
                            group_of = NULL) {
  walk(file.path(out_dir, c("b_reports", "c_data")), dir.create,
    recursive = TRUE, showWarnings = FALSE
  )
  read_sheet <- \(sheet) readxl::read_excel(paths[["contrasts"]], sheet)
  set_tests <- readRDS(paths[["set_tests"]])
  proteins <- set_tests$protein_results
  phenotype <- readr::read_csv(paths[["phenotype"]], show_col_types = FALSE)
  n_proteins <- n_distinct(proteins$uniprot_id)

  counts <- panel_protein_counts(proteins, groups, n_proteins, validity)
  pca <- panel_pca_pair(
    read_sheet("pca"), read_sheet("pca_variance")$variance_explained, read_sheet("permanova"),
    phenotype, groups, group_of, n_proteins
  )
  pathways <- panel_pathway_grid(set_tests, readRDS(paths[["gene_sets"]]), groups)
  fry_rows <- read_sheet("fry_concordance")
  barcodes <- map(
    set_names(c("Training", "Acute")),
    \(window) signature_barcode(proteins, fry_rows, window, groups = groups)
  )

  counted <- summarise(counts$data, fdr = sum(fdr), .by = contrast) |>
    inner_join(read_sheet("contrast_summary"), by = "contrast")
  barcode_numbers <- map(barcodes, "labels") |> list_rbind()
  fry_numbers <- fry_rows |>
    filter(startsWith(set, "nominal_"), sub("_.*", "", set_from) == sub("_.*", "", ranked_on)) |>
    transmute(
      window = sub("_.*", "", set_from), chosen = sub(".*_", "", set_from),
      ranked_in = sub(".*_", "", ranked_on), direction = sub(".*_", "", set), fdr_source = fdr
    )
  traced <- inner_join(barcode_numbers, fry_numbers,
    by = c("window", "chosen", "ranked_in", "direction")
  )
  stopifnot(
    nrow(counted) == 9, all(counted$fdr == counted$n_fdr_05),
    all(map_int(pca$tests, nrow) == 1),
    all(count(pathways$data, selected_in, set_id)$n == pathways$n_columns),
    nrow(traced) == 8, all(traced$fdr == traced$fdr_source)
  )

  rows <- list(
    # free() keeps a and c from aligning to the long labels in panel d.
    (patchwork::free(panel_a) | counts$plot | patchwork::free(pca$plot)) +
      patchwork::plot_layout(widths = c(1.1, 1.45, 2)),
    (pathways$plot | barcodes$Training$plot | barcodes$Acute$plot) +
      patchwork::plot_layout(widths = c(1.45, 1, 1))
  )
  figure <- patchwork::wrap_plots(rows, ncol = 1, heights = c(66, 98)) +
    # The phenotype PCA shares tag c with the proteome PCA.
    patchwork::plot_annotation(tag_levels = list(c("a", "b", "c", "", "d", "e", "f"))) &
    figure_tags
  save_figure(figure, file.path(out_dir, "b_reports"), stem, width = 300, height = 175)

  sheets <- c(a_sheets, list(
    protein_counts = counts$data,
    pca_proteome = pca$proteome,
    pca_phenotype = pca$phenotype,
    permanova = read_sheet("permanova"),
    pathways = pathways$data,
    barcode = select(
      barcode_numbers, window, chosen, ranked_in, direction, n, same_sign_share, fdr
    )
  ))
  holds <- c(
    a_holds,
    "Panel b: per contrast and direction, proteins at Pi < 0.05, at p < 0.05 or Pi, at FDR < 0.05.",
    "Panel c, left: PC1 and PC2 per sample (protein-contrasts stage, pca).",
    "Panel c, right: PC1 and PC2 per biopsy from the six scaled size and strength traits.",
    "Panel c: PERMANOVA, proteome (time, group) and phenotype (time), protein-contrasts stage.",
    "Panel d: each selected set in all nine contrasts, with NES and fgsea FDR.",
    "Panels e, f: each group's p < 0.05 sets in the other group's ranking: n, % same way, fry FDR."
  )
  write_workbook(
    sheets, holds, file.path(out_dir, "c_data", paste0(stem, ".xlsx")), inputs, paths,
    note = note
  )
}

# Panel a for a cluster arc: each person on the first two components of what they were clustered
# on, coloured by group, with HR/LR as shape; people not clustered are grey. The subtitle gives the
# K the gap statistic selects, the K used, and each group's size and bootstrap Jaccard.
panel_clustering <- function(space, variance, choice, stability, groups, clustered_on) {
  d <- space |>
    mutate(group = factor(group, groups), arm = factor(arm, c("HR", "LR")))
  jaccard <- stability |>
    filter(k == choice$k_used) |>
    inner_join(distinct(filter(d, !is.na(group)), cluster, group), by = "cluster") |>
    arrange(group)
  sizes <- count(filter(d, !is.na(group)), group)
  axis <- \(i) sprintf("PC%d (%.0f%%)", i, 100 * variance$variance_explained[i])
  plot <- ggplot(d, aes(PC1, PC2)) +
    geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey85") +
    geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey85") +
    geom_point(aes(colour = group, shape = arm), size = 2) +
    scale_colour_manual(values = group_colours(groups), na.value = "grey80", name = NULL) +
    scale_shape_manual(values = arm_shapes, name = NULL) +
    labs(
      title = "Who is in each group",
      subtitle = sprintf(
        "PAM on %s; gap selects K = %d, K = %d used\n%s",
        clustered_on, choice$k_gap, choice$k_used,
        paste(sprintf(
          "%s n = %d, Jaccard %.2f", jaccard$group,
          sizes$n[match(jaccard$group, sizes$group)], jaccard$jaccard
        ), collapse = "; ")
      ),
      x = axis(1), y = axis(2)
    ) +
    figure_theme +
    theme(legend.position = "bottom", legend.box = "vertical", legend.spacing.y = unit(0, "mm"))
  list(plot = plot, data = select(d, subject, arm, group, cluster, PC1, PC2), jaccard = jaccard)
}

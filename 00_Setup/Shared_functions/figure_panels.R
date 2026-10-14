# Panel builders shared by the figures in every Figures folder. Each takes its data as
# arguments and returns a ggplot (or a list holding one and the numbers it plots); colours and
# theme come from figure_theme.R, which must be sourced first.

# PCA scores with a shaded 95% t-ellipse per time point and a test summary in the top-left corner.
pca_panel <- function(scores, variance, title, subtitle, stats, shapes = arm_shapes) {
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
    ggplot2::scale_shape_manual(values = shapes, name = NULL) +
    # Headroom above the points for the test summary.
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.32))) +
    ggplot2::labs(x = axis(1), y = axis(2), title = title, subtitle = subtitle) +
    figure_theme
}

# One window's signature rank: each group's sets for that window (p < 0.05 or Pi < 0.05) placed in
# the other group's ranking for the same window, drawn as limma's barcode (Wu & Smyth 2012).
# Members are ticks, the up-set above the axis and the down-set below, with a tricube moving
# average of membership along the rank. The two groups are different people, so the comparison
# shares no group mean and fry's result (02_Protein_Contrasts, fry_concordance) is a clean test.
signature_barcode <- function(proteins, fry_rows, window, set_level = "nominal",
                              span = 0.12, floor = 0.34, groups = c("HR", "LR")) {
  contrasts <- paste0(window, "_", groups)
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
      uniprot_id, ranked_in = ifelse(chosen == groups[1], groups[2], groups[1])
    ) |>
    dplyr::inner_join(rankings, by = c("ranked_in", "uniprot_id")) |>
    dplyr::summarise(ranks = list(rank), n = dplyr::n(), .by = c(chosen, ranked_in, direction)) |>
    dplyr::mutate(
      lane = factor(lane_name(chosen, ranked_in), lane_name(rev(groups), groups)),
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
        "rank in the other group's %s contrast, up to down (n = %d)", tolower(window), n_ranked
      ),
      y = NULL, title = paste(window, "signature rank"),
      subtitle = sprintf(
        "fry: each group's %s %s sets in the other group's ranking", level_text, tolower(window)
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

# Keeps each set unless it shares at least cut of its measured proteins with a set already kept;
# ids must be ordered strongest first.
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

# PDF only, through cairo, which embeds fonts and carries the Pi glyph the base device drops.
save_figure <- function(plot, dir, stem, width, height) {
  ggplot2::ggsave(
    file.path(dir, paste0(stem, ".pdf")), plot,
    width = width, height = height, units = "mm", device = grDevices::cairo_pdf
  )
}


# The builders below serve F2 and F3 in every grouping. Each returns list(plot, data, source): the
# plot, the numbers it draws, and one row naming where they come from. They expect ggplot2, dplyr,
# tidyr and purrr attached, and pick rows by the rule each states.
panel_source <- \(panel, shows, file, sheet, rule) {
  tibble::tibble(panel, shows, file, sheet, rule)
}

# Strip boxes for facets named after something other than a contrast family, in one family's tint,
# or in grey when family is NA.
family_strips <- function(family, n, axis = "x") {
  tone <- \(name, grey) dplyr::coalesce(contrast_tone(name)[family], grey)
  boxes <- ggh4x::elem_list_rect(
    fill = rep(tone("fill", "grey92"), n), colour = rep(tone("line", "grey60"), n),
    linewidth = 0.6
  )
  text <- ggh4x::elem_list_text(colour = "black", face = "bold")
  if (axis == "x") {
    ggh4x::strip_themed(background_x = boxes, text_x = text)
  } else {
    ggh4x::strip_themed(background_y = boxes, text_y = text)
  }
}

# A module's name for figures: its top over-represented set at FDR < alpha_fdr
# (07_Module_Characterisation), or a statement that it has none.
module_identity <- function(labels, width = 34) {
  identity <- dplyr::if_else(
    labels$top_set_fdr < alpha_fdr, pathway_label(labels$top_set, width), "no set at FDR < 0.05"
  )
  rlang::set_names(paste0(identity, " (", labels$module, ")"), labels$module)
}

# What each module is: its size in its own WGCNA colour, named by its top set, with its three
# three top members by kME. Rows run largest first, the order panel_module_contrasts() uses.
panel_module_identity <- function(labels, file) {
  d <- labels |>
    mutate(
      module = factor(module, rev(module)),
      top_members = map_chr(
        strsplit(top_members, ", "), \(h) paste(utils::head(h, 3), collapse = ", ")
      )
    )
  plot <- ggplot(d, aes(n_proteins, module)) +
    geom_col(aes(fill = module), width = 0.6, colour = "grey30", linewidth = 0.15, alpha = 0.85) +
    geom_text(
      aes(label = top_members),
      hjust = -0.08, size = 1.7, colour = "grey35", fontface = "italic"
    ) +
    scale_fill_identity() +
    scale_y_discrete(labels = module_identity(labels)) +
    scale_x_continuous(expand = expansion(mult = c(0, 1.1)), breaks = c(0, 100, 200)) +
    facet_grid2(. ~ "Module", strip = family_strips(NA, 1)) +
    labs(
      title = "Co-expression programmes",
      subtitle = sprintf("%d WGCNA modules: size, top set, top members", nrow(labels)),
      x = "proteins", y = NULL
    ) +
    figure_theme +
    theme(panel.grid.major.y = element_blank(), axis.text.y = element_text(size = 5.5))
  list(
    plot = plot,
    data = select(d, module, n_proteins, top_set, top_set_fdr, top_members, string_ratio),
    source = panel_source(
      "a", "module size, top set and top members", file, "labels",
      "every module; top set if FDR < 0.05"
    )
  )
}

# Every module on every contrast: tile fill is the eigengene logFC; the dot is limma's evidence on
# the eigengene, filled at FDR and open at p only; the square marks fry p < alpha_nominal on the
# module's proteins. Rows follow module_levels, so the panel sits beside panel_module_identity().
panel_module_contrasts <- function(eigengene_tests, fry_tests, module_levels, groups, file) {
  d <- eigengene_tests |>
    select(contrast, module, logFC, p, fdr) |>
    left_join(
      select(fry_tests, contrast, module, p_fry = p, fdr_fry = fdr),
      by = c("contrast", "module")
    ) |>
    mutate(
      family = contrast_panel(contrast),
      column = factor(contrast_short(contrast), unique(contrast_short(contrast))),
      module = factor(module, rev(module_levels)),
      eigengene = evidence_level(p, fdr)
    )
  marks <- bind_rows(
    filter(d, eigengene != evidence_levels[3]) |>
      mutate(mark = paste("eigengene", sub(" only", "", eigengene))),
    filter(d, p_fry < alpha_nominal) |> mutate(mark = sprintf("fry p < %s", alpha_nominal))
  )
  mark_shapes <- c(16, 1, 0) |>
    set_names(c(
      sprintf("eigengene FDR < %s", alpha_fdr), sprintf("eigengene p < %s", alpha_nominal),
      sprintf("fry p < %s", alpha_nominal)
    ))
  limit <- ceiling(max(abs(d$logFC)) * 20) / 20
  plot <- ggplot(d, aes(column, module)) +
    geom_tile(aes(fill = logFC), colour = "white", linewidth = 0.3) +
    geom_point(aes(shape = mark, size = mark), marks, stroke = 0.35) +
    scale_fill_distiller(
      palette = "PiYG", direction = 1, limits = c(-limit, limit), name = "eigengene\nlogFC"
    ) +
    scale_shape_manual(values = mark_shapes, name = NULL) +
    scale_size_manual(values = set_names(c(1, 1, 2.4), names(mark_shapes)), name = NULL) +
    facet_grid2(
      . ~ family,
      scales = "free_x", space = "free_x", strip = contrast_strips(x = levels(d$family))
    ) +
    labs(
      title = "Programmes across the nine contrasts",
      subtitle = sprintf(
        "%s minus %s for between-group and interaction columns", groups[1], groups[2]
      ),
      x = NULL, y = NULL
    ) +
    figure_theme +
    theme(
      axis.text.y = element_blank(), panel.grid.major = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.spacing.x = unit(0.8, "mm"), legend.position = "right",
      legend.box = "vertical", strip.text.x = element_text(size = 6, face = "bold")
    ) +
    guides(fill = guide_colourbar(order = 1, barheight = unit(18, "mm")))
  list(
    plot = plot, data = select(d, -eigengene),
    source = panel_source(
      "b", "eigengene logFC and evidence, fry evidence", file, "eigengene_tests, fry_tests",
      "every module and contrast"
    )
  )
}

# Each person's eigengene at T1, T2 and T3 for one module, one facet per group, with the group mean;
# the subtitle gives the module's eigengene test in each interaction contrast.
panel_module_people <- function(eigengenes, meta, module, tests, identity, groups, file, rule) {
  d <- tibble::tibble(sample_id = colnames(eigengenes), value = eigengenes[module, ]) |>
    inner_join(select(meta, sample_id, subject, arm, timepoint), by = "sample_id") |>
    mutate(arm = factor(arm, groups))
  means <- summarise(d, value = mean(value), .by = c(arm, timepoint))
  stats <- tests |>
    filter(.data$module == .env$module, grepl("Interaction", contrast)) |>
    with(paste(sprintf("%s logFC %.2f, %s", contrast_label(contrast), logFC, label_p(p)),
      collapse = "; "
    ))
  plot <- ggplot(d, aes(timepoint, value, colour = arm)) +
    geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey70") +
    geom_line(aes(group = subject), alpha = 0.5, linewidth = 0.25) +
    geom_line(aes(group = arm), data = means, linewidth = 1) +
    geom_point(data = means, size = 1.4) +
    facet_wrap(~arm, nrow = 1) +
    scale_colour_manual(values = group_colours(groups), guide = "none") +
    scale_x_discrete(expand = expansion(add = 0.2)) +
    labs(
      title = sprintf("%s, person by person", identity[[module]]),
      subtitle = paste0(stats, "\nThin lines: people; thick: group mean"),
      x = NULL, y = "eigengene"
    ) +
    figure_theme +
    theme(strip.text = element_text(face = "bold", size = 6.5))
  list(
    plot = plot, data = select(d, subject, arm, timepoint, value) |> mutate(module = module),
    source = panel_source("d", "eigengene per person", file, "eigengenes", rule)
  )
}

# Gene sets called per contrast by fgsea (FDR), fry (FDR) and fry (nominal p), up and down. fgsea
# ranks a set against the other proteins; fry asks whether it moved away from zero across people,
# so fry needs more people behind a shift. The dashed line is the fry count expected by chance.
panel_set_calls <- function(set_tests, contrast_order, groups, file) {
  calls <- c(
    sprintf("fgsea\nFDR < %s", alpha_fdr), sprintf("fry\nFDR < %s", alpha_fdr),
    sprintf("fry\np < %s", alpha_nominal)
  )
  tests <- filter(set_tests, method %in% c("fgsea", "fry"))
  n_tested <- n_distinct(tests$set_id)
  counts <- bind_rows(
    filter(tests, method == "fgsea", fdr < alpha_fdr) |> mutate(call = calls[1]),
    filter(tests, method == "fry", fdr < alpha_fdr) |> mutate(call = calls[2]),
    filter(tests, method == "fry", p < alpha_nominal) |> mutate(call = calls[3])
  ) |>
    count(call, contrast, direction) |>
    complete(
      call = calls, contrast = contrast_order, direction = c("Up", "Down"), fill = list(n = 0L)
    ) |>
    mutate(
      call = factor(call, calls), family = contrast_panel(contrast),
      row = factor(contrast_short(contrast), rev(unique(contrast_short(contrast_order)))),
      direction = factor(direction, c("Up", "Down"))
    )
  totals <- summarise(counts, n = sum(n), .by = c(call, family, row, contrast))
  plot <- ggplot(counts, aes(n, row)) +
    geom_col(aes(fill = direction), width = 0.7, position = position_stack(reverse = TRUE)) +
    geom_text(aes(label = n), totals, hjust = -0.15, size = 1.7) +
    geom_vline(
      aes(xintercept = x),
      tibble::tibble(call = factor(calls[3], calls), x = n_tested * alpha_nominal),
      linetype = "dashed", linewidth = 0.3
    ) +
    facet_grid2(
      family ~ call,
      scales = "free_y", space = "free_y", switch = "y",
      strip = contrast_strips(y = levels(counts$family))
    ) +
    scale_fill_manual(values = direction_colours, name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.35)), n.breaks = 3) +
    force_panelsizes(cols = c(1, 0.3, 1)) +
    labs(
      title = "Gene sets per contrast, fgsea and fry",
      subtitle = sprintf(
        paste(
          "fgsea ranks a set against the other proteins; fry tests its change across people.",
          "\n%s sets; dashed: %s fry calls expected by chance"
        ),
        scales::comma(n_tested), round(n_tested * alpha_nominal)
      ),
      x = NULL, y = NULL
    ) +
    figure_theme +
    theme(
      panel.grid.major = element_blank(), axis.text.x = element_blank(),
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 90, size = 5.5), legend.position = "bottom",
      panel.spacing.y = unit(1, "mm"), panel.spacing.x = unit(3, "mm")
    )
  list(
    plot = plot, data = totals,
    source = panel_source(
      "c", "sets called per contrast", file, "set_tests", "every set and contrast"
    )
  )
}

# The five between-group classification tasks of
# 03_Group_comparisons/a_script/05_classification.qmd under the contrast strips: the
# levels sit under Between groups, the changes under their window.
task_panel <- function(task) {
  family <- dplyr::case_when(
    startsWith(task, "Training_change") ~ "Training",
    startsWith(task, "Acute_change") ~ "Acute",
    .default = "Between groups"
  )
  factor(family, panel_families)
}
task_short <- function(task) {
  short <- c(
    Baseline = "T1", Trained = "T2", Acute = "T3", Training_change = "T2 - T1",
    Acute_change = "T3 - T2"
  )
  unname(short[sub("_[^_]*$", "", task)])
}

# The best single feature per task and level (smallest Wilcoxon p, then AUC furthest from 0.5)
# with its 95% interval: pROC's stratified bootstrap where the stage drew the feature, DeLong's
# otherwise. The right-hand count is features at p < alpha_nominal against the number chance gives.
panel_auc_top <- function(auc_by_level, counts_by_level, groups, file) {
  best <- imap(auc_by_level, \(auc, level) {
    auc |>
      filter(!is.na(p)) |>
      arrange(p, desc(abs(auc - 0.5))) |>
      slice_head(n = 1, by = task) |>
      mutate(level = level)
  }) |>
    list_rbind() |>
    mutate(
      interval = if_else(is.na(boot_lower), "DeLong", "bootstrap"),
      lower = coalesce(boot_lower, auc_lower), upper = coalesce(boot_upper, auc_upper)
    )
  counts <- imap(counts_by_level, \(counts, level) {
    summarise(counts, n_nominal = sum(n_nominal), n_tested = sum(n_tested), .by = task) |>
      mutate(level = level)
  }) |>
    list_rbind()
  d <- inner_join(best, counts, by = c("level", "task")) |>
    mutate(
      level = factor(level, names(auc_by_level)), family = task_panel(task),
      row = factor(task_short(task), rev(task_short(unique(task)))),
      name = if_else(level == "sets", pathway_label(sub(".*[|]", "", feature), 18), label),
      tally = sprintf("%s  %d / %d", name, n_nominal, round(n_tested * alpha_nominal))
    )
  plot <- ggplot(d, aes(auc, row)) +
    geom_vline(xintercept = 0.5, linetype = "dashed", linewidth = 0.3, colour = "grey50") +
    geom_errorbar(
      aes(xmin = lower, xmax = upper, linetype = interval),
      orientation = "y", width = 0.3, linewidth = 0.35
    ) +
    geom_point(aes(colour = auc > 0.5), size = 1.5) +
    geom_text(aes(x = 1.05, label = tally), hjust = 0, size = 1.7, colour = "grey20") +
    facet_grid2(
      family ~ level,
      scales = "free_y", space = "free_y", switch = "y",
      strip = contrast_strips(y = levels(d$family)),
      labeller = labeller(level = stringr::str_to_sentence)
    ) +
    scale_colour_manual(
      values = set_names(unname(arm_colours), c("TRUE", "FALSE")), breaks = c("TRUE", "FALSE"),
      labels = paste("higher in", groups), name = NULL
    ) +
    scale_linetype_manual(values = c(bootstrap = "solid", DeLong = "dotted"), name = "95% CI") +
    scale_x_continuous(breaks = c(0, 0.5, 1)) +
    coord_cartesian(xlim = c(0, 1.85), clip = "off") +
    labs(
      title = sprintf("Can one feature tell %s from %s?", groups[1], groups[2]),
      subtitle = sprintf(
        paste(
          "Best feature per task (smallest Wilcoxon p), AUC with 95%% CI;",
          "right: that feature, then features at p < %s / expected by chance"
        ),
        alpha_nominal
      ),
      x = "AUC", y = NULL
    ) +
    figure_theme +
    theme(
      strip.placement = "outside", strip.text.y.left = element_text(angle = 90, size = 5.5),
      strip.text.x = element_text(face = "bold"), panel.grid.major.y = element_blank(),
      legend.position = "bottom", panel.spacing.x = unit(4, "mm")
    )
  list(
    plot = plot,
    data = select(
      d, level, task, feature, label, auc, auc_lower, auc_upper, boot_lower, boot_upper, p, fdr,
      n_positive, n_negative, n_nominal, n_tested
    ),
    source = panel_source(
      "a", "best single feature per task", file, "auc, counts",
      "smallest p per task and level, then AUC furthest from 0.5"
    )
  )
}

# The group-by-trait difference for every protein in one window and trait (t of the first group's
# slope minus the second's) against the continuous growth-by-strength interaction (t) for the same
# protein. The group split is a split of growth, so the two codings ask nearly one question.
panel_group_vs_continuous <- function(group_difference, interaction, window, outcome, groups,
                                      file, n_label = 8) {
  d <- group_difference |>
    filter(level == "proteins", .data$window == .env$window, .data$trait == .env$outcome) |>
    select(feature, label, t_group = t, p_group = p, fdr_group = fdr) |>
    inner_join(
      filter(interaction, .data$window == .env$window, strength_test == outcome) |>
        select(feature, t_continuous = t, p_continuous = p, fdr_continuous = fdr),
      by = "feature"
    ) |>
    mutate(hit = evidence_level(p_group, fdr_group))
  rho <- stats::cor(d$t_group, d$t_continuous, method = "spearman")
  labelled <- slice_min(filter(d, fdr_group < alpha_fdr), p_group, n = n_label)
  plot <- ggplot(d, aes(t_group, t_continuous)) +
    geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey60") +
    geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey60") +
    geom_point(data = filter(d, hit == evidence_levels[3]), size = 0.4, colour = "grey80") +
    geom_point(
      aes(colour = hit),
      data = filter(d, hit != evidence_levels[3]), size = 0.9
    ) +
    ggrepel::geom_text_repel(
      aes(label = label), labelled,
      size = 1.8, seed = 1, min.segment.length = 0, segment.size = 0.2, max.overlaps = Inf
    ) +
    scale_colour_manual(
      values = set_names(
        contrast_tone("ink")[c("Interaction", "Between groups")], evidence_levels[1:2]
      ),
      name = sprintf("%s - %s", groups[1], groups[2])
    ) +
    labs(
      title = "Group difference against the continuous interaction",
      subtitle = sprintf(
        "%s, %s window, %d proteins; Spearman across proteins %.2f (descriptive)",
        trait_label(outcome), window, nrow(d), rho
      ),
      x = sprintf("t, %s minus %s slope (group x trait)", groups[1], groups[2]),
      y = "t, growth x strength (continuous)"
    ) +
    figure_theme +
    theme(legend.position = "bottom")
  list(
    plot = plot, data = d,
    source = panel_source(
      "b", "group difference against continuous interaction", file,
      "group_difference; interaction", sprintf("proteins, %s, %s", window, outcome)
    )
  )
}

# Each protein's within-group Spearman rho with one trait in the first group against the second,
# with the group-difference hits marked. Within-group rho is descriptive at six to eight people.
panel_within_rho <- function(within_group, group_difference, window, outcome, groups, file) {
  hits <- group_difference |>
    filter(level == "proteins", .data$window == .env$window, .data$trait == .env$outcome) |>
    transmute(feature, hit = evidence_level(p, fdr))
  d <- within_group |>
    filter(level == "proteins", .data$window == .env$window, .data$trait == .env$outcome) |>
    select(feature, label, group, rho) |>
    pivot_wider(names_from = group, values_from = rho) |>
    rename(first = all_of(groups[1]), second = all_of(groups[2])) |>
    inner_join(hits, by = "feature")
  n_crossed <- sum(d$hit == evidence_levels[1] & sign(d$first) != sign(d$second))
  n_hits <- sum(d$hit == evidence_levels[1])
  plot <- ggplot(d, aes(first, second)) +
    geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey60") +
    geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey60") +
    geom_abline(linetype = "dashed", linewidth = 0.2, colour = "grey60") +
    geom_point(data = filter(d, hit != evidence_levels[1]), size = 0.4, colour = "grey80") +
    geom_point(
      data = filter(d, hit == evidence_levels[1]), size = 1.1,
      colour = contrast_tone("ink")[["Interaction"]]
    ) +
    ggrepel::geom_text_repel(
      aes(label = label), filter(d, hit == evidence_levels[1]),
      size = 1.8, seed = 1, min.segment.length = 0, segment.size = 0.2, max.overlaps = Inf
    ) +
    coord_equal(xlim = c(-1, 1), ylim = c(-1, 1)) +
    labs(
      title = "The same proteins, within each group",
      subtitle = sprintf(
        "Spearman rho with %s, %s window; %d of %d group-difference hits change sign",
        trait_label(outcome), window, n_crossed, n_hits
      ),
      x = sprintf("rho in %s", groups[1]), y = sprintf("rho in %s", groups[2])
    ) +
    figure_theme
  list(
    plot = plot, data = d |> rename_with(\(x) paste0("rho_", groups), c(first, second)),
    source = panel_source(
      "c", "within-group rho", file, "within_group; group_difference",
      sprintf("proteins, %s, %s; marked: group-difference FDR < %s", window, outcome, alpha_fdr)
    )
  )
}

# One protein's change against one trait, person by person, a least-squares line per group, and
# each group's Spearman rho from the within-group sheet.
panel_group_slope <- function(row, values, outcome, group_of, within, groups, file) {
  rho <- within |>
    filter(feature == row$feature, window == row$window, trait == row$trait) |>
    with(set_names(sprintf("%s: rho %.2f (n = %d)", group, rho, n_subjects), group))
  d <- tibble::tibble(subject = names(values), value = values, trait = outcome[names(values)]) |>
    filter(!is.na(trait), !is.na(value), subject %in% names(group_of)) |>
    mutate(arm = factor(group_of[subject], groups))
  plot <- ggplot(d, aes(value, trait, colour = arm)) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.6) +
    geom_point(aes(shape = arm), size = 1.4) +
    scale_colour_manual(values = group_colours(groups), labels = rho[groups], name = NULL) +
    scale_shape_manual(
      values = set_names(unname(arm_shapes), groups), labels = rho[groups], name = NULL
    ) +
    labs(
      title = sprintf("%s, person by person", row$label),
      subtitle = sprintf(
        "Group difference in slope: t = %.1f, %s, %s", row$t, label_p(row$p), label_fdr(row$fdr)
      ),
      x = sprintf("%s, %s change (log2)", row$label, row$window), y = trait_label(row$trait)
    ) +
    figure_theme +
    theme(legend.position = "bottom", legend.direction = "vertical")
  list(
    plot = plot, data = d,
    source = panel_source(
      "d", "lead group-difference protein per person", file, "group_difference, within_group",
      "smallest group-difference p among proteins"
    )
  )
}

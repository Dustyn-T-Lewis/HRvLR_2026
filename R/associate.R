# The continuous arc's association step, shared by the protein, set and module levels. Each level's
# script reads its own matrix and labels and passes them here. No group label enters: every subject
# is pooled and each trait is a number.
#
# Sample model: every T1 and T2 biopsy carries its own trait values. Per trait, limma fits
# abundance ~ timepoint + between + within with subject blocked, where between is the subject's
# mean trait and within the biopsy's deviation from it. within asks whether a feature rises in the
# biopsies where that person's trait rose; between asks whether people with more of the trait
# carry more of the feature.
#
# Change score: one value per subject and pairing. Pairings: training change, T1 level and acute
# change against the % trait changes, composite and volume ("baseline" is the T1 level against
# later change); T1 and T2 levels against the trait levels at the same biopsy; and the T3 level
# against the T2 trait levels, since T3 follows the same training (it also carries the bout's
# acute effect). Per feature: limma on one value per subject (p), Spearman (p_rank) and biweight
# midcorrelation (bicor).
#
# Two predictor: feature ~ growth + strength gain. Interaction: feature ~ growth x strength gain,
# the continuous form of asking whether the strength link depends on growth.
#
# collection_of, when given, makes each collection its own screen: BH runs within collection.
associate_step <- function(stage, inputs, paths, values, meta, feature_label, feature_noun,
                           feature_unit, value_label, row_groups, row_order,
                           collection_of = NULL, panel_name = \(label, collection) label,
                           observed_share = NULL, save_change_score = FALSE) {
  out <- file.path(stage, "c_data")
  figures <- file.path(stage, "b_reports")
  step <- tolower(basename(stage))
  workbook <- paste0(step, ".xlsx")
  for (path in c(out, figures)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  by_collection <- if (is.null(collection_of)) character(0) else "collection"
  labelled_heat <- \(...) heat(..., labels = feature_label)
  labelled <- function(d) {
    if (length(by_collection)) d <- mutate(d, collection = collection_of[feature], .after = feature)
    d <- mutate(d, label = feature_label[feature], .after = all_of(c("feature", by_collection)))
    if (is.null(observed_share)) {
      return(d)
    }
    mutate(d, observed_share = observed_share[feature], .after = label)
  }
  corrected <- function(d, ...) {
    mutate(d, fdr = p.adjust(p, "BH"), qvalue = storey_q(p), .by = c(..., all_of(by_collection)))
  }

  phenotype <- readr::read_csv(paths[["phenotype"]], show_col_types = FALSE)
  phenotype_list <- phenotype_tables(phenotype)
  traits <- phenotype_list$traits
  biopsy <- phenotype_list$biopsy
  outcomes <- phenotype_list$outcomes
  pairing_tables <- phenotype_list$pairings
  pairing_outcomes <- map(pairing_tables, \(tt) setdiff(names(tt), "subject"))
  outcome_of <- \(name, subjects) pairing_value(pairing_tables, "training", name, subjects)

  samples <- meta |>
    filter(timepoint %in% c("T1", "T2")) |>
    select(sample_id, subject, timepoint) |>
    left_join(select(biopsy, -arm), by = c("subject", "timepoint"))
  model <- map(set_names(traits), \(trait) {
    s <- samples |>
      filter(!is.na(.data[[trait]])) |>
      mutate(
        between = ave(.data[[trait]], subject),
        within = .data[[trait]] - between
      )
    design <- model.matrix(~ timepoint + between + within, s)
    y <- values[, s$sample_id]
    correlation <- duplicateCorrelation(y, design, block = s$subject)$consensus.correlation
    fit <- lmFit(y, design, block = s$subject, correlation = correlation) |>
      eBayes(robust = TRUE)
    map(c("within", "between"), \(term) {
      topTable(fit, coef = term, number = Inf, sort.by = "none") |>
        as_tibble(rownames = "feature") |>
        transmute(trait, term, feature, n_samples = nrow(s), effect = logFC, t, p = P.Value)
    }) |>
      list_rbind()
  }) |>
    list_rbind() |>
    labelled() |>
    corrected(trait, term)

  windows <- list(
    training = subject_change(values, meta, "T1", "T2"),
    baseline = subject_values(values, meta, "T1"),
    acute = subject_change(values, meta, "T2", "T3"),
    baseline_level = subject_values(values, meta, "T1"),
    post_level = subject_values(values, meta, "T2"),
    post_bout_level = subject_values(values, meta, "T3")
  )
  change_windows <- c("training", "baseline", "acute")
  linear <- \(x, outcome) subject_fit(x, data.frame(trait = outcome), ~trait, "trait")
  continuous <- function(x, outcome) {
    fit <- linear(x, outcome) |> select(feature, n_subjects, effect, t, p)
    ok <- !is.na(outcome)
    bc <- WGCNA::bicorAndPvalue(t(x[fit$feature, ok, drop = FALSE]), outcome[ok])
    fit |>
      left_join(select(spearman_by_row(x, outcome), feature, rho, p_rank = p), by = "feature") |>
      mutate(bicor = bc$bicor[feature, 1], p_bicor = bc$p[feature, 1])
  }
  change_score <- imap(windows, \(x, window) {
    map(set_names(pairing_outcomes[[window]]), \(name) {
      continuous(x, pairing_value(pairing_tables, window, name, colnames(x)))
    }) |>
      list_rbind(names_to = "outcome") |>
      mutate(window = window, .before = 1)
  }) |>
    list_rbind() |>
    labelled() |>
    corrected(window, outcome)
  strength_tests <- c("pct_leg_press_1rm", "pct_leg_ext_1rm")
  two_predictor <- imap(windows[change_windows], \(x, window) {
    map(set_names(strength_tests), \(test) {
      predictors <- data.frame(
        growth = outcome_of("comp_hypertrophy", colnames(x)),
        strength = outcome_of(test, colnames(x))
      )
      subject_fit(x, predictors, ~ growth + strength, c("growth", "strength"))
    }) |>
      list_rbind(names_to = "strength_test") |>
      mutate(window = window, .before = 1)
  }) |>
    list_rbind() |>
    mutate(
      outcome = if_else(coefficient == "growth", "comp_hypertrophy", strength_test),
      .after = coefficient
    ) |>
    labelled() |>
    corrected(window, strength_test, coefficient)

  # Interaction: does the link to strength gain depend on how much the person grew?
  interaction <- imap(windows[change_windows], \(x, window) {
    map(set_names(strength_tests), \(test) {
      predictors <- data.frame(
        growth = outcome_of("comp_hypertrophy", colnames(x)),
        strength = outcome_of(test, colnames(x))
      )
      subject_fit(x, predictors, ~ growth * strength, "growth:strength")
    }) |>
      list_rbind(names_to = "strength_test") |>
      mutate(window = window, .before = 1)
  }) |>
    list_rbind() |>
    select(-coefficient) |>
    labelled() |>
    corrected(window, strength_test)

  tally <- function(d, ...) {
    summarise(d,
      n_tested = sum(!is.na(p)), n_nominal = sum(p < 0.05, na.rm = TRUE),
      n_fdr_05 = sum(fdr < 0.05, na.rm = TRUE), .by = c(..., all_of(by_collection))
    )
  }
  counts <- bind_rows(
    tally(model, term, trait) |>
      transmute(
        analysis = paste("model:", term), comparison = trait, pick(all_of(by_collection)),
        n_tested, n_nominal, n_fdr_05
      ),
    tally(change_score, window, outcome) |>
      transmute(
        analysis = paste("change score:", window), comparison = outcome,
        pick(all_of(by_collection)), n_tested, n_nominal, n_fdr_05
      ),
    tally(two_predictor, window, strength_test, coefficient) |>
      transmute(
        analysis = paste("two predictor:", window),
        comparison = paste(coefficient, strength_test, sep = ", "), pick(all_of(by_collection)),
        n_tested, n_nominal, n_fdr_05
      ),
    tally(interaction, window, strength_test) |>
      transmute(
        analysis = paste("interaction:", window), comparison = strength_test,
        pick(all_of(by_collection)), n_tested, n_nominal, n_fdr_05
      )
  )

  # Units come from the phenotype column names.
  unit_cols <- grep("_pre_", names(phenotype), value = TRUE)
  unit_of <- set_names(sub(".*_pre_", "", unit_cols), sub("_pre_.*", "", unit_cols))
  outcome_unit <- c(
    comp_hypertrophy = "score", volume_load_total_kg = "kg",
    set_names(rep("%", length(traits)), paste0("pct_", traits)), unit_of[traits]
  )
  window_axis <- c(
    training = "change, T1 to T2", baseline = "T1 level", acute = "change, T2 to T3",
    baseline_level = "T1 level", post_level = "T2 level", post_bout_level = "T3 level"
  )
  bh_scope <- if (length(by_collection)) " and collection" else ""
  name_of <- \(d) panel_name(d$label, d[["collection"]])

  # Overviews: sample-model t by trait, then change-score rho by outcome, one page per window.
  t_limit <- cap(model$t)
  overviews <- imap(row_groups, \(rows, group) {
    term_grid <- \(term, value) {
      as_grid(filter(model, term == .env$term), rows, traits, "trait", value)
    }
    model_maps <- labelled_heat(
      term_grid("within", \(d) d$t), "Within-person t", c(-t_limit, 0, t_limit),
      "t"
    ) +
      labelled_heat(
        term_grid("between", \(d) d$t), "Between-person t", c(-t_limit, 0, t_limit), NULL
      ) +
      labelled_heat(term_grid("within", \(d) status_of(d$p, d$fdr)), "Within: p and BH q") +
      labelled_heat(term_grid("between", \(d) status_of(d$p, d$fdr)), "Between: p and BH q",
        row_names = TRUE
      )
    model_page <- heatmap_page(model_maps, sprintf("Sample model overview: %s", group), paste(
      sprintf("All %s %s, one row each.", format(length(rows), big.mark = ","), feature_noun),
      "Moderated t of the within and between terms, abundance ~ timepoint + between + within,",
      sprintf("subject blocked. Colours stop at the 99th percentile of |t|, %.2g.", t_limit),
      row_order, paste0("Grey: not tested. Table: c_data/", workbook, ", model.")
    ))
    window_pages <- map(names(windows), \(w) {
      d <- filter(change_score, window == w)
      maps <- labelled_heat(
        as_grid(d, rows, pairing_outcomes[[w]], "outcome", \(d) d$rho),
        "Spearman rho, both arms pooled", c(-1, 0, 1), "rho"
      ) +
        labelled_heat(
          as_grid(d, rows, pairing_outcomes[[w]], "outcome", \(d) status_of(d$p, d$fdr)),
          "limma p and BH q",
          row_names = TRUE
        )
      heatmap_page(maps, sprintf("Change score overview, %s window: %s", w, group), paste(
        sprintf("All %s %s, one row each.", format(length(rows), big.mark = ","), feature_noun),
        sprintf(
          "Spearman rho between the %s (%s) and each outcome; status from limma.", value_label,
          window_axis[[w]]
        ),
        row_order, paste0("Grey: not tested. Table: c_data/", workbook, ", change_score.")
      ))
    })
    c(list(model_page), window_pages)
  }) |>
    list_flatten()

  # Grids: every p < 0.05 hit, strongest first (p, then |t| or |rho|, then feature ID).
  trait_samples <- map(set_names(traits), \(trait) filter(samples, !is.na(.data[[trait]])))
  seen_in <- \(feature, trait) !is.na(values[feature, trait_samples[[trait]]$sample_id])
  model_hits <- model |>
    filter(p < 0.05) |>
    mutate(
      n_biopsies = map2_int(feature, trait, \(f, tr) sum(seen_in(f, tr))),
      n_subjects = map2_int(feature, trait, \(f, tr) {
        n_distinct(trait_samples[[tr]]$subject[seen_in(f, tr)])
      }),
      name = name_of(pick(everything())),
      panel = sprintf(
        "%s\n%s %.3g  p %.2g  q %.2g\n%d biopsies, %d subjects", name, term, effect, p, fdr,
        n_biopsies, n_subjects
      )
    ) |>
    arrange(match(trait, traits), desc(term), p, desc(abs(t)), feature) |>
    paginate(c("trait", "term"))
  model_scope <- if (length(by_collection)) "trait, term and collection," else "trait and term,"
  model_page <- function(hits) {
    trait <- hits$trait[1]
    term <- hits$term[1]
    info <- counts |>
      filter(analysis == paste("model:", term), comparison == trait) |>
      summarise(across(c(n_tested, n_nominal, n_fdr_05), sum))
    points <- hits |>
      select(feature, panel) |>
      cross_join(select(samples, sample_id, subject, timepoint, trait_value = all_of(trait))) |>
      mutate(
        value = values[cbind(feature, sample_id)],
        panel = factor(panel, slots(hits$panel))
      ) |>
      filter(!is.na(value), !is.na(trait_value))
    ggplot(points, aes(trait_value, value)) +
      geom_line(
        aes(group = subject),
        data = \(d) filter(d, n() == 2, .by = c(panel, subject)), alpha = 0.6, colour = "grey40"
      ) +
      geom_point(aes(shape = timepoint), size = 1.5) +
      facet_wrap(~panel, nrow = 3, ncol = 4, scales = "free", drop = FALSE) +
      scale_shape_manual(values = c(T1 = 1, T2 = 16), name = NULL) +
      labs(
        title = sprintf("%s, %s-person term", trait, term),
        subtitle = counts_line(info, feature_noun, hits),
        x = sprintf("%s at the biopsy (%s)", trait, unit_of[[trait]]), y = value_label,
        caption = stringr::str_wrap(width = 160, paste(
          "limma, abundance ~ timepoint + between + within, subject blocked, on the T1 and T2",
          "biopsies with the trait measured. Lines join each subject's two biopsies: the within",
          "term reads their slopes, the between term the spread of subjects. Header: effect per",
          "trait unit, p, BH q within", model_scope,
          "and the biopsies and subjects where the feature was observed.",
          paste0("Table: c_data/", workbook, ", model.")
        ))
      ) +
      page_theme
  }

  change_hits <- change_score |>
    filter(p < 0.05) |>
    mutate(
      name = name_of(pick(everything())),
      panel = sprintf(
        "%s vs %s (%s)\nrho %+.2f  p %.2g  q %.2g  n %d", name, outcome, outcome_unit[outcome],
        rho, p, fdr, n_subjects
      )
    ) |>
    arrange(match(window, names(windows)), p, desc(abs(rho)), feature, outcome) |>
    paginate("window")
  page_captions <- c(
    "change score" = paste(
      "One point per subject, every subject pooled. The line is a least-squares fit, drawn as a",
      "guide. Header: Spearman rho, and limma p and BH q within window, outcome%s, with the",
      paste0("subjects used. Table: c_data/", workbook, ", change_score.")
    ),
    "two predictor" = paste(
      "One point per subject. Header: the growth or strength coefficient from limma, feature ~",
      "comp_hypertrophy + strength gain, each adjusted for the other, with p, BH q within",
      "window, strength test, term%s, and the subjects used. The x axis shows the term's own",
      "outcome; the drawing does not hold the other predictor fixed.",
      paste0("Table: c_data/", workbook, ", two_predictor.")
    )
  )
  page_captions[] <- sprintf(page_captions, bh_scope)
  change_page <- function(hits, kind = "change score") {
    w <- hits$window[1]
    info <- counts |>
      filter(analysis == paste0(kind, ": ", w)) |>
      summarise(across(c(n_tested, n_nominal, n_fdr_05), sum))
    points <- hits |>
      select(feature, outcome, panel) |>
      cross_join(tibble(subject = colnames(windows[[w]]))) |>
      mutate(
        value = windows[[w]][cbind(feature, subject)],
        outcome_value = map2_dbl(outcome, subject, \(o, s) pairing_value(pairing_tables, w, o, s)),
        panel = factor(panel, slots(hits$panel))
      ) |>
      filter(!is.na(value), !is.na(outcome_value))
    ggplot(points, aes(value, outcome_value)) +
      # A straight line needs two points, not the default 80.
      geom_smooth(
        method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.5, n = 2, colour = "grey40"
      ) +
      geom_point(size = 1.3) +
      facet_wrap(~panel, nrow = 3, ncol = 4, scales = "free", drop = FALSE) +
      labs(
        title = sprintf("%s, %s window", stringr::str_to_sentence(kind), w),
        subtitle = counts_line(info, "feature-outcome pairs", hits),
        x = sprintf("%s, %s", value_label, window_axis[[w]]), y = "outcome (unit in header)",
        caption = stringr::str_wrap(page_captions[[kind]], 160)
      ) +
      page_theme
  }
  model_pages <- map(split(model_hits, model_hits$pdf_page), model_page)
  change_pages <- map(split(change_hits, change_hits$pdf_page), change_page)
  slope_hits <- two_predictor |>
    filter(p < 0.05) |>
    mutate(
      kind = "two predictor",
      header = sprintf(
        "%s term, %s model %+.3g", coefficient, sub("pct_", "", strength_test), effect
      )
    ) |>
    mutate(
      name = name_of(pick(everything())),
      panel = sprintf(
        "%s vs %s (%s)\n%s  p %.2g  q %.2g  n %d", name, outcome, outcome_unit[outcome], header,
        p, fdr, n_subjects
      )
    ) |>
    arrange(match(window, names(windows)), p, desc(abs(t)), feature, outcome) |>
    paginate(c("kind", "window"))
  slope_pages <- map(split(slope_hits, slope_hits$pdf_page), \(hits) {
    change_page(hits, hits$kind[1])
  })

  plot_index <- bind_rows(
    tibble(
      analysis = "overview", label = map_chr(overviews, "title"), pdf_page = seq_along(overviews)
    ),
    transmute(model_hits,
      analysis = "model", trait, term, pick(all_of(by_collection)), feature, label,
      pdf_page = pdf_page + length(overviews), position
    ),
    transmute(change_hits,
      analysis = "change score", window, outcome, pick(all_of(by_collection)), feature, label,
      pdf_page = pdf_page + length(overviews) + length(model_pages), position
    ),
    transmute(slope_hits,
      analysis = kind, window, outcome, strength_test, coefficient, pick(all_of(by_collection)),
      feature, label,
      pdf_page = pdf_page + length(overviews) + length(model_pages) + length(change_pages),
      position
    )
  )
  stopifnot(
    nrow(model_hits) == sum(model$p < 0.05, na.rm = TRUE),
    nrow(change_hits) == sum(change_score$p < 0.05),
    !anyDuplicated(select(filter(plot_index, !is.na(feature)), -label, -pdf_page, -position)),
    nrow(slope_hits) == sum(two_predictor$p < 0.05),
    max(plot_index$pdf_page) ==
      length(overviews) + length(model_pages) + length(change_pages) + length(slope_pages)
  )

  if (save_change_score) {
    saveRDS(list(change_score = change_score), file.path(out, "associate.rds"), compress = "xz")
  }

  # Points as one Dingbats glyph each rather than four curves: the PDF shrinks by about 60%.
  pdf(file.path(figures, paste0(step, "_figures.pdf")),
    width = 11, height = 8.5, useDingbats = TRUE
  )
  walk(overviews, \(page) page$draw())
  walk(c(model_pages, change_pages, slope_pages), print)
  invisible(dev.off())

  write_workbook(
    list(
      counts = counts, model = arrange(model, trait, term, p),
      change_score = arrange(change_score, window, outcome, p),
      two_predictor = arrange(two_predictor, window, strength_test, coefficient, p),
      interaction = arrange(interaction, window, strength_test, p),
      plot_index = plot_index
    ),
    c(
      "Tested, p < 0.05 and BH < 0.05 counts per term and trait, and per window and outcome.",
      sprintf(
        "Sample model: effect per trait unit, t, p, BH and Storey q per %s, trait and term.",
        feature_unit
      ),
      sprintf(paste(
        "Change score, every subject pooled: limma effect, t, p, BH, Storey q; Spearman rho and",
        "p_rank; bicor and its p; per %s, window and outcome."
      ), feature_unit),
      sprintf(
        "Two predictor: growth and strength terms, t, p, BH and Storey q per %s, window and test.",
        feature_unit
      ),
      sprintf(paste(
        "Interaction, feature ~ growth x strength gain: effect, t, p, BH and Storey q per %s,",
        "window and test."
      ), feature_unit),
      "PDF page of each overview, and page and panel position of every p < 0.05 row."
    ),
    file.path(out, workbook),
    inputs, paths
  )
}

# The 04_Associate step, shared by the protein, set and module levels. Each level's script reads its
# own matrix and labels and passes them here.
#
# Sample model: every T1 and T2 biopsy carries its own trait values. Per trait, limma fits
# abundance ~ timepoint + between + within with subject blocked, where between is the subject's
# mean trait and within the biopsy's deviation from it. within asks whether a feature rises in the
# biopsies where that person's trait rose; between asks whether people with more of the trait
# carry more of the feature. Arm stays out: it was defined from one of these traits.
#
# Change score: one value per subject and pairing, every subject pooled, each trait continuous.
# Pairings: training change, T1 level and acute change against the % trait changes, composite and
# volume ("baseline" is the T1 level against later change), and T1 and T2 levels against the trait
# levels at the same biopsy. Per feature: limma on one value per subject (p, the primary test),
# Spearman (p_rank), biweight midcorrelation (bicor), and the largest limma p over leave-one-out
# refits; robust = all four below 0.05. Each screen's count and strongest hit are judged against
# 200 shuffles of whole phenotype rows (sheet chance), with the empirical FDR = shuffle median /
# observed count. Within-arm Spearman p is exact.
#
# Interaction: feature change ~ growth x strength gain, with the leave-one-out maximum p and the
# person whose removal gives it.
# Flips: DGCA's Fisher z difference between HR and LR correlations, p from 200 arm-label
# permutations.
#
# collection_of, when given, makes each collection its own screen: BH runs within collection.
associate_step <- function(stage, inputs, paths, values, meta, feature_label, feature_noun,
                           feature_unit, value_label, row_groups, row_order,
                           collection_of = NULL, panel_name = \(label, collection) label,
                           observed_share = NULL,
                           all_arm_rows = FALSE, save_change_score = FALSE) {
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
  arm_of <- set_names(phenotype$arm, phenotype$subject)
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
    post_level = subject_values(values, meta, "T2")
  )
  change_windows <- c("training", "baseline", "acute")
  linear <- \(x, outcome) subject_fit(x, data.frame(trait = outcome), ~trait, "trait")
  continuous <- function(x, outcome) {
    fit <- linear(x, outcome) |> select(feature, n_subjects, effect, t, p)
    ok <- !is.na(outcome)
    bc <- WGCNA::bicorAndPvalue(t(x[fit$feature, ok, drop = FALSE]), outcome[ok])
    loo <- map(colnames(x)[ok], \(s) {
      keep <- colnames(x) != s
      linear(x[, keep, drop = FALSE], outcome[keep]) |> select(feature, p)
    }) |>
      list_rbind() |>
      summarise(loo_max_p = max(p), .by = feature)
    fit |>
      left_join(select(spearman_by_row(x, outcome), feature, rho, p_rank = p), by = "feature") |>
      mutate(bicor = bc$bicor[feature, 1], p_bicor = bc$p[feature, 1]) |>
      left_join(loo, by = "feature") |>
      mutate(robust = p < 0.05 & p_rank < 0.05 & p_bicor < 0.05 & loo_max_p < 0.05)
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
  by_arm <- map(set_names(c("HR", "LR")), \(arm) {
    imap(windows[change_windows], \(x, window) {
      kept <- x[, arm_of[colnames(x)] == arm, drop = FALSE]
      map(set_names(outcomes), \(name) {
        spearman_by_row(kept, outcome_of(name, colnames(kept)), exact = TRUE)
      }) |>
        list_rbind(names_to = "outcome") |>
        mutate(window = window, .before = 1)
    }) |>
      list_rbind()
  }) |>
    list_rbind(names_to = "arm") |>
    filter(!is.na(rho))

  # Chance line per screen: whole phenotype rows shuffled across subjects, so traits keep their
  # correlations and only the link to the proteome breaks.
  null_screens <- parallel::mclapply(1:200, \(i) {
    shuffled <- shuffle_pairings(pairing_tables, i)
    imap(windows, \(x, window) {
      map(set_names(pairing_outcomes[[window]]), \(name) {
        d <- linear(x, pairing_value(shuffled, window, name, colnames(x)))
        if (length(by_collection)) d <- mutate(d, collection = collection_of[feature])
        summarise(d, nominal = sum(p < 0.05), min_p = min(p), .by = any_of("collection"))
      }) |>
        list_rbind(names_to = "outcome") |>
        mutate(window = window)
    }) |>
      list_rbind()
  }, mc.cores = 8) |>
    list_rbind()
  chance <- change_score |>
    summarise(
      tested = n(), nominal = sum(p < 0.05), bh_05 = sum(fdr < 0.05),
      robust = sum(robust, na.rm = TRUE), min_p = min(p),
      .by = c(window, outcome, all_of(by_collection))
    ) |>
    left_join(
      summarise(null_screens,
        null_median = median(nominal),
        null_95 = quantile(nominal, 0.95), null_counts = list(nominal), null_min = list(min_p),
        .by = c(window, outcome, all_of(by_collection))
      ),
      by = c("window", "outcome", by_collection)
    ) |>
    mutate(
      empirical_fdr = pmin(1, null_median / nominal),
      perm_p_count = map2_dbl(null_counts, nominal, \(v, o) (sum(v >= o) + 1) / (length(v) + 1)),
      perm_p_strongest = map2_dbl(null_min, min_p, \(v, o) (sum(v <= o) + 1) / (length(v) + 1))
    ) |>
    select(-null_counts, -null_min)

  # comp_hypertrophy defines the arms, so a slope difference on it re-measures the split; it
  # stays in the pooled correlation and the two-predictor model. Each arm's slope needs five
  # subjects with both the feature and the outcome observed.
  arm_difference <- imap(windows[change_windows], \(x, window) {
    map(set_names(setdiff(outcomes, "comp_hypertrophy")), \(name) {
      predictors <- data.frame(
        arm = factor(arm_of[colnames(x)], c("LR", "HR")), outcome = outcome_of(name, colnames(x))
      )
      seen <- !is.na(x) & rep(!is.na(predictors$outcome), each = nrow(x))
      per_arm <- sapply(c("HR", "LR"), \(a) rowSums(seen[, predictors$arm %in% a, drop = FALSE]))
      subject_fit(
        x[apply(per_arm, 1, min) >= 5, , drop = FALSE], predictors, ~ arm * outcome,
        "armHR:outcome"
      )
    }) |>
      list_rbind(names_to = "outcome") |>
      mutate(window = window, .before = 1)
  }) |>
    list_rbind() |>
    select(-coefficient) |>
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
      fit <- subject_fit(x, predictors, ~ growth * strength, "growth:strength")
      ok <- complete.cases(predictors)
      loo <- map(colnames(x)[ok], \(s) {
        keep <- colnames(x) != s
        subject_fit(
          x[, keep, drop = FALSE], predictors[keep, ], ~ growth * strength,
          "growth:strength"
        ) |>
          transmute(feature, p, left_out = s)
      }) |>
        list_rbind() |>
        # The person whose removal weakens the result most is the one it leans on.
        summarise(loo_max_p = max(p), loo_worst_subject = left_out[which.max(p)], .by = feature)
      left_join(fit, loo, by = "feature")
    }) |>
      list_rbind(names_to = "strength_test") |>
      mutate(window = window, .before = 1)
  }) |>
    list_rbind() |>
    select(-coefficient) |>
    labelled() |>
    corrected(window, strength_test)

  # Flips: does a feature's correlation with a trait differ between HR and LR? DGCA on the change
  # windows and the T1 and T2 levels; each arm needs five subjects with the trait measured.
  flips <- imap(windows[c("training", "acute", "baseline_level", "post_level")], \(x, window) {
    map(setdiff(pairing_outcomes[[window]], "comp_hypertrophy"), \(name) {
      outcome <- pairing_value(pairing_tables, window, name, colnames(x))
      ok <- !is.na(outcome)
      arm <- factor(arm_of[colnames(x)][ok], c("HR", "LR"))
      if (min(table(arm)) < 5) {
        return(NULL)
      }
      trait_row <- matrix(outcome[ok], nrow = 1, dimnames = list(name, NULL))
      design <- DGCA::makeDesign(arm)
      arm_data <- x[, ok, drop = FALSE]
      observed <- DGCA::ddcorAll(arm_data, design, c("HR", "LR"),
        inputMatB = trait_row, corrType = "spearman", adjust = "none", sigThresh = 1,
        nPairs = "all"
      ) |>
        suppressMessages()
      set.seed(1)
      null <- DGCA::getDCorPerm(arm_data, design, c("HR", "LR"),
        inputMatB = trait_row, nPerms = 200, corrType = "spearman"
      ) |>
        suppressMessages() |>
        as.vector() |>
        abs()
      null <- sort(null[!is.na(null)])
      # Pooled over features as DGCA does; left.open counts an infinite z (a perfect correlation
      # in a small arm) only against null values that are also infinite.
      empirical <- \(z) {
        (1 + length(null) - findInterval(abs(z), null, left.open = TRUE)) / (1 + length(null))
      }
      as_tibble(observed) |>
        transmute(
          feature = as.character(Gene1), outcome = name, hr_rho = HR_cor,
          lr_rho = LR_cor, z_diff = if_else(is.finite(zScoreDiff), zScoreDiff, NA_real_),
          p = empirical(zScoreDiff), class = Classes
        )
    }) |>
      list_rbind() |>
      mutate(window = window, .before = 1)
  }) |>
    list_rbind() |>
    labelled() |>
    corrected(window, outcome)

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
    tally(arm_difference, window, outcome) |>
      transmute(
        analysis = paste("arm difference:", window), comparison = outcome,
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
      ),
    tally(flips, window, outcome) |>
      transmute(
        analysis = paste("flips:", window), comparison = outcome, pick(all_of(by_collection)),
        n_tested, n_nominal, n_fdr_05
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
    baseline_level = "T1 level", post_level = "T2 level"
  )
  arm_colours <- c(HR = "#2166AC", LR = "#B2182B")
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
        value = values[cbind(feature, sample_id)], arm = arm_of[subject],
        panel = factor(panel, slots(hits$panel))
      ) |>
      filter(!is.na(value), !is.na(trait_value))
    ggplot(points, aes(trait_value, value, colour = arm)) +
      geom_line(
        aes(group = subject),
        data = \(d) filter(d, n() == 2, .by = c(panel, subject)), alpha = 0.6
      ) +
      geom_point(aes(shape = timepoint), size = 1.5) +
      facet_wrap(~panel, nrow = 3, ncol = 4, scales = "free", drop = FALSE) +
      scale_colour_manual(values = arm_colours, name = NULL) +
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

  arm_text <- by_arm |>
    mutate(text = sprintf("%-3s n=%d rho=%+.2f p=%.3f", arm, n_subjects, rho, p)) |>
    summarise(caption = paste(text, collapse = "\n"), .by = c(window, outcome, feature))
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
      "One point per subject. Lines are a least-squares fit per arm, drawn as a guide. Header:",
      "Spearman rho, and limma p and BH q pooled over both arms within window, outcome%s, with",
      "the subjects used.",
      "Corner: each arm's rho and exact p. Arm is the median split of",
      "comp_hypertrophy, so pooled correlations with it and with fibre area change partly",
      paste0("carry the group split. Table: c_data/", workbook, ", change_score and"),
      "change_score_by_arm."
    ),
    "arm difference" = paste(
      "One point per subject. Lines are a least-squares fit per arm. Header: HR's slope minus",
      "LR's per outcome unit from limma, feature ~ arm * outcome, with p, BH q within window,",
      "outcome%s, and the subjects used. Corner: each arm's Spearman rho and exact p.",
      paste0("Table: c_data/", workbook, ", arm_difference.")
    ),
    "two predictor" = paste(
      "One point per subject. Header: the growth or strength coefficient from limma, feature ~",
      "comp_hypertrophy + strength gain, each adjusted for the other, with p, BH q within",
      "window, strength test, term%s, and the subjects used. The x axis shows the term's own",
      "outcome; the drawing does not hold the other predictor fixed. Corner: each arm's",
      paste0("Spearman rho and exact p. Table: c_data/", workbook, ", two_predictor.")
    )
  )
  page_captions[] <- sprintf(page_captions, bh_scope)
  change_page <- function(hits, kind = "change score") {
    w <- hits$window[1]
    info <- counts |>
      filter(analysis == paste0(kind, ": ", w)) |>
      summarise(across(c(n_tested, n_nominal, n_fdr_05), sum))
    points <- hits |>
      left_join(arm_text, by = c("window", "outcome", "feature")) |>
      select(feature, outcome, panel, caption) |>
      cross_join(tibble(subject = colnames(windows[[w]]))) |>
      mutate(
        value = windows[[w]][cbind(feature, subject)],
        outcome_value = map2_dbl(outcome, subject, \(o, s) pairing_value(pairing_tables, w, o, s)),
        arm = arm_of[subject], caption = coalesce(caption, ""),
        panel = factor(panel, slots(hits$panel))
      ) |>
      filter(!is.na(value), !is.na(outcome_value))
    ggplot(points, aes(value, outcome_value, colour = arm)) +
      # A straight line needs two points, not the default 80.
      geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.5, n = 2) +
      geom_point(size = 1.3) +
      geom_text(
        data = distinct(points, panel, caption), aes(-Inf, Inf, label = caption),
        inherit.aes = FALSE, hjust = -0.03, vjust = 1.1, size = 1.9, family = "mono"
      ) +
      facet_wrap(~panel, nrow = 3, ncol = 4, scales = "free", drop = FALSE) +
      scale_colour_manual(values = arm_colours, name = NULL) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.45))) +
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
  slope_hits <- bind_rows(
    arm_difference |>
      filter(p < 0.05) |>
      mutate(kind = "arm difference", header = sprintf("slope HR-LR %+.3g", effect)),
    two_predictor |>
      filter(p < 0.05) |>
      mutate(
        kind = "two predictor",
        header = sprintf(
          "%s term, %s model %+.3g", coefficient, sub("pct_", "", strength_test), effect
        )
      )
  ) |>
    mutate(
      name = name_of(pick(everything())),
      panel = sprintf(
        "%s vs %s (%s)\n%s  p %.2g  q %.2g  n %d", name, outcome, outcome_unit[outcome], header,
        p, fdr, n_subjects
      )
    ) |>
    arrange(
      match(kind, c("arm difference", "two predictor")), match(window, names(windows)), p,
      desc(abs(t)), feature, outcome
    ) |>
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
    nrow(slope_hits) == sum(arm_difference$p < 0.05) + sum(two_predictor$p < 0.05),
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
      chance = chance,
      change_score_by_arm = if (all_arm_rows) by_arm else filter(by_arm, p < 0.05),
      arm_difference = arrange(arm_difference, window, outcome, p),
      two_predictor = arrange(two_predictor, window, strength_test, coefficient, p),
      interaction = arrange(interaction, window, strength_test, p),
      flips = arrange(flips, window, outcome, p),
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
        "p_rank; bicor; leave-one-out maximum p; robust; per %s, window and outcome."
      ), feature_unit),
      paste(
        "Per window and outcome: tested, p < 0.05, BH < 0.05, robust; the median and 95th",
        "percentile over 200 phenotype shuffles, empirical FDR, and permutation p for the count",
        "and for the strongest hit."
      ),
      if (all_arm_rows) {
        sprintf("Change score within each arm, exact p, every %s.", feature_unit)
      } else {
        "Change score within each arm, exact p, rows with p < 0.05."
      },
      sprintf(
        "Arm difference: HR minus LR slope, t, p, BH and Storey q per %s, window and outcome.",
        feature_unit
      ),
      sprintf(
        "Two predictor: growth and strength terms, t, p, BH and Storey q per %s, window and test.",
        feature_unit
      ),
      sprintf(paste(
        "Interaction, feature ~ growth x strength gain: effect, t, p, BH, Storey q, the",
        "leave-one-out maximum p and the person left out to give it, per %s, window and test."
      ), feature_unit),
      sprintf(paste(
        "Flips, DGCA: HR and LR Spearman rho, Fisher z difference (blank where one arm's rho is",
        "+1 or -1), p from 200 arm-label permutations, BH and Storey q per %s, window and outcome."
      ), feature_unit),
      "PDF page of each overview, and page and panel position of every p < 0.05 row."
    ),
    file.path(out, workbook),
    inputs, paths
  )
}

# Run once per proteome clustering (T1, T2, T3, training change, acute change) on that
# clustering's groups C1 and C2: a copy of 02_HRvLR_groups/02_DEP. Each input writes to its own
# folder; stack_inputs() then joins them into c_data/all_inputs.xlsx with a validity tag.
# The nine contrasts on limpa's quantification: limma with precision weights from limpa's standard
# errors, subject blocked, blood as a covariate; each logFC with its 95% confidence interval, BH
# within each contrast, null checks, and fry on one arm's responders in the other arm.

suppressPackageStartupMessages({
  library(here)
  library(limma)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(grid)
  library(ComplexHeatmap)
  library(writexl)
})
source(here("R", "helpers.R"))
source(here("03_Cluster_groups", "Proteome", "R", "inputs.R"))

run_input <- function(input) {
  stage <- here("03_Cluster_groups", "Proteome", "03_DEP", input)
  out <- file.path(stage, "c_data")
  figures <- file.path(stage, "b_reports")
  for (path in c(out, figures)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  inputs <- c(
    design = "03_Cluster_groups/Proteome/02_Design/{input}/c_data/design.rds",
    proteins = "03_Cluster_groups/Proteome/01_Cluster/c_data/inputs/{input}/proteins.rds",
    corrected = "03_Cluster_groups/Proteome/01_Cluster/c_data/inputs/{input}/DAList_corrected.rds",
    phenotype = "00_Input/phenotype.csv"
  )
  inputs <- sub("{input}", input, inputs, fixed = TRUE)
  paths <- map_chr(inputs, here)
  stopifnot(file.exists(paths))
  d <- readRDS(paths[["design"]])
  dal <- d$dal
  contrast_names <- colnames(dal$design$contrast_matrix)
  stopifnot(identical(rownames(dal$design$design_matrix), colnames(dal$data)))

  proteins <- readRDS(paths[["proteins"]])$proteins
  design <- dal$design$design_matrix
  contrast_matrix <- dal$design$contrast_matrix
  meta <- dal$metadata
  stopifnot(identical(colnames(proteins$E), rownames(design)))
  eb <- fit_limpa(proteins, design, contrast_matrix, meta$subject)

  # The same nine contrasts without the blood term: the sensitivity check for the blood covariate.
  blood_free <- design[, colnames(design) != "blood"]
  raw_eb <- fit_limpa(proteins, blood_free, contrast_matrix[colnames(blood_free), ], meta$subject)
  unadjusted <- tibble(
    uniprot_id = rep(rownames(proteins$E), length(contrast_names)),
    contrast = rep(contrast_names, each = nrow(proteins$E)),
    logFC = c(raw_eb$coefficients[, contrast_names]), p = c(raw_eb$p.value[, contrast_names])
  ) |>
    pivot_wider(
      names_from = contrast, values_from = c(logFC, p), names_glue = "{contrast}_{.value}"
    )

  # Effect the primary contrast could detect at 80% power, two-sided 0.05.
  se <- eb$stdev.unscaled[, "Training_Interaction"] * sqrt(eb$s2.post)
  detectable <- se * (qt(0.975, eb$df.total) + qt(0.80, eb$df.total))

  annotation <- as_tibble(dal$annotation)[, c("uniprot_id", "gene", "protein", "description")]
  stopifnot(identical(rownames(eb$coefficients), annotation$uniprot_id))
  # limpa gives every protein a value in every sample, so every contrast tests every protein.
  # pi_score (Xiao et al. 2014, inverted) only selects sets for the fry check below.
  results <- map(set_names(contrast_names), \(contrast) {
    topTable(eb, coef = contrast, number = Inf, sort.by = "none", confint = TRUE) |>
      as_tibble() |>
      select(logFC, ci_lower = CI.L, ci_upper = CI.R, AveExpr, t, P.Value, adj.P.Val) |>
      bind_cols(annotation, x = _)
  }) |>
    list_rbind(names_to = "contrast") |>
    rename(p = P.Value, fdr = adj.P.Val) |>
    mutate(qvalue = storey_q(p), .by = contrast) |>
    mutate(
      pi_score = p^abs(logFC),
      sig_pi = case_when(
        pi_score < 0.05 & logFC > 0 ~ 1L,
        pi_score < 0.05 & logFC < 0 ~ -1L,
        TRUE ~ 0L
      )
    )
  n_obs <- set_names(rowSums(proteins$other$n.observations > 0), rownames(proteins$E))

  contrast_summary <- results |>
    summarise(
      n_tested = sum(!is.na(p)),
      n_nominal = sum(p < 0.05, na.rm = TRUE),
      n_fdr_05 = sum(fdr < 0.05, na.rm = TRUE),
      n_fdr_10 = sum(fdr < 0.10, na.rm = TRUE),
      min_fdr = min(fdr, na.rm = TRUE),
      .by = contrast
    ) |>
    left_join(d$roles, by = "contrast")

  # propTrueNull() saturates at 1, so it cannot tell a flat null from a conservative one. Median |t|
  # (0.674 under a calibrated null) and the share of p below 0.2 can.
  null_calibration <- results |>
    filter(!is.na(p)) |>
    summarise(
      prop_true_null = round(propTrueNull(p), 3),
      median_abs_t = round(median(abs(t)), 3),
      p_below_0.2 = round(mean(p < 0.2), 3), .by = contrast
    )

  # fry on one arm's responders, tested in another contrast: the other arm's same window, and each
  # arm's acute bout for the training sets. Only the other-arm pairs share no cell mean; a training
  # set tested on its own arm's acute contrast shares T2, so a protein high at T2 by chance reads as
  # up then down. fry reads limpa's complete matrix with the same design and block.
  abundance <- proteins$E
  stopifnot(
    identical(colnames(abundance), rownames(design)), identical(meta$sample_id, rownames(design)),
    !anyNA(abundance)
  )
  fry_correlation <- duplicateCorrelation(abundance, design, block = meta$subject)
  fry_correlation <- fry_correlation$consensus.correlation
  member_rows <- function(cn, criterion, direction) {
    r <- filter(results, contrast == cn)
    sign_ok <- if (direction == "up") r$logFC > 0 else r$logFC < 0
    hit <- switch(criterion,
      pi = r$sig_pi == if (direction == "up") 1L else -1L,
      bh = !is.na(r$fdr) & r$fdr < 0.05 & sign_ok,
      nominal = !is.na(r$p) & r$p < 0.05 & sign_ok
    )
    which(rownames(abundance) %in% r$uniprot_id[coalesce(hit, FALSE)])
  }
  pairings <- tribble(
    ~set_from, ~ranked_on,
    "Training_C1", "Training_C2",
    "Training_C2", "Training_C1",
    "Acute_C1", "Acute_C2",
    "Acute_C2", "Acute_C1",
    "Training_C1", "Acute_C1",
    "Training_C1", "Acute_C2",
    "Training_C2", "Acute_C2",
    "Training_C2", "Acute_C1"
  )
  selections <- expand_grid(criterion = c("pi", "bh", "nominal"), direction = c("up", "down"))
  fry_concordance <- pmap(pairings, \(set_from, ranked_on) {
    sets <- pmap(selections, \(criterion, direction) member_rows(set_from, criterion, direction)) |>
      set_names(paste(selections$criterion, selections$direction, sep = "_"))
    sets <- keep(sets, \(rows) length(rows) >= 5)
    # fry() drops its FDR column when handed a single set.
    if (length(sets) < 2) {
      return(NULL)
    }
    # The effect sizes behind the test: each set's mean log2 FC in the arm that chose it (inflated
    # by the choice) and in the other arm (different people, so free of it), and the share of
    # members whose sign carries over.
    log_fc <- \(cn, rows) filter(results, contrast == cn)$logFC[rows]
    carried <- imap(sets, \(rows, set) {
      own <- log_fc(set_from, rows)
      other <- log_fc(ranked_on, rows)
      tibble(
        set,
        mean_logfc_own = mean(own), mean_logfc_other = mean(other),
        same_sign_share = mean(sign(own) == sign(other))
      )
    }) |>
      list_rbind()
    fry(abundance,
      index = sets, design = design, contrast = contrast_matrix[, ranked_on],
      block = meta$subject, correlation = fry_correlation
    ) |>
      as_tibble(rownames = "set") |>
      transmute(
        set_from, ranked_on, set,
        n_proteins = NGenes, direction = Direction, p = PValue, fdr = FDR,
        p_mixed = PValue.Mixed, fdr_mixed = FDR.Mixed
      ) |>
      left_join(carried, by = "set")
  }) |>
    list_rbind()

  p_histogram <- results |>
    filter(!is.na(p)) |>
    mutate(contrast = factor(contrast, levels = contrast_names)) |>
    ggplot(aes(p)) +
    geom_histogram(binwidth = 0.05, boundary = 0, fill = "grey35") +
    geom_hline(
      aes(yintercept = n_tested / 20),
      mutate(contrast_summary, contrast = factor(contrast, levels = contrast_names)),
      linetype = 2, colour = "firebrick"
    ) +
    facet_wrap(~contrast) +
    labs(
      title = "p-value histogram per contrast", x = "p", y = "proteins",
      caption = paste(
        "Dashed line: a flat null (tested / 20). A spike at zero on a flat background is signal;",
        "a slope up toward 1 is a conservative test. Table: c_data/02_contrasts.xlsx, DEP_matrix."
      )
    ) +
    theme_minimal(base_size = 10)

  # One volcano per contrast. Nominal proteins in colour; BH < 0.05 proteins labelled.
  volcano <- results |>
    filter(!is.na(p)) |>
    mutate(
      contrast = factor(contrast, levels = contrast_names),
      call = case_when(fdr < 0.05 ~ "BH < 0.05", p < 0.05 ~ "p < 0.05", TRUE ~ "not nominal")
    ) |>
    ggplot(aes(logFC, -log10(p), colour = call)) +
    geom_point(size = 0.6, alpha = 0.7) +
    geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "grey50") +
    ggrepel::geom_text_repel(
      aes(label = gene),
      data = \(x) filter(x, fdr < 0.05), size = 2.2, show.legend = FALSE
    ) +
    facet_wrap(~contrast) +
    scale_colour_manual(
      values = c(`BH < 0.05` = "#B2182B", `p < 0.05` = "#2166AC", `not nominal` = "grey75"),
      name = NULL
    ) +
    labs(
      title = "Volcano per contrast",
      subtitle = sprintf(
        "%s proteins at p < 0.05 per contrast; %d at BH < 0.05 across all nine.",
        paste(range(contrast_summary$n_nominal), collapse = " to "), sum(contrast_summary$n_fdr_05)
      ),
      x = "log2 fold change", y = "-log10 p",
      caption = paste(
        "limma on limpa's quantification, subject blocked, blood as a covariate. Dashed line:",
        "p = 0.05. BH within contrast. Table: c_data/02_contrasts.xlsx, DEP_matrix."
      )
    ) +
    theme_minimal(base_size = 9) +
    theme(legend.position = "top", plot.caption = element_text(hjust = 0))

  # Later steps read the fit, its design and limpa's matrix from here.
  saveRDS(
    list(
      eBayes_fit = eb, design = dal$design, metadata = meta, data = proteins$E,
      annotation = dal$annotation, proteins = proteins
    ),
    file.path(out, "fit.rds"),
    compress = "xz"
  )
  sample_values <- proteins$E
  column_meta <- meta
  feature_label <- set_names(annotation$gene, annotation$uniprot_id)
  feature_noun <- "proteins"
  contrast_results <- transmute(results, feature = uniprot_id, contrast, effect = logFC, p, fdr)
  # Display order: average-linkage clustering on 1 - Pearson r of limpa's matrix.
  row_groups <- list(
    proteins = rownames(abundance)[hclust(as.dist(1 - cor(t(abundance))), "average")$order]
  )

  # Overviews: every feature across the 45 biopsies, then across the nine contrasts. Biopsy traits
  # sit on the T1 and T2 columns; T3 has none. comp_hypertrophy is one value per subject.
  arm_colours <- c(C1 = "#2166AC", C2 = "#B2182B")
  trait_table <- readr::read_csv(paths[["phenotype"]], show_col_types = FALSE)
  biopsy_traits <- trait_table |>
    select(subject, matches("_(pre|post)_")) |>
    pivot_longer(-subject, names_to = c("trait", "when"), names_pattern = "(.*)_(pre|post)_.*") |>
    pivot_wider(names_from = trait) |>
    transmute(
      subject,
      timepoint = if_else(when == "pre", "T1", "T2"),
      fcsa_mixed, mcsa, leg_press_1rm
    )
  columns <- as_tibble(column_meta) |>
    transmute(sample_id, subject, arm = as.character(arm), timepoint = as.character(timepoint)) |>
    mutate(group = factor(paste(arm, timepoint, sep = "_"), paste(
      rep(c("C1", "C2"), each = 3), c("T1", "T2", "T3"),
      sep = "_"
    ))) |>
    left_join(biopsy_traits, by = c("subject", "timepoint")) |>
    left_join(select(trait_table, subject, comp_hypertrophy), by = "subject") |>
    arrange(group, subject)
  stopifnot(!anyDuplicated(columns$sample_id), setequal(columns$sample_id, colnames(sample_values)))
  grey_ramp <- \(x) circlize::colorRamp2(range(x, na.rm = TRUE), c("#F7F7F7", "#252525"))
  tracks <- list(
    "fCSA mixed (um2)" = columns$fcsa_mixed,
    "mCSA (cm2)" = columns$mcsa, "leg press 1RM (kg)" = columns$leg_press_1rm,
    "comp_hypertrophy, per subject" = columns$comp_hypertrophy
  )
  sample_annotation <- HeatmapAnnotation(
    df = data.frame(arm = columns$arm, timepoint = columns$timepoint, tracks, check.names = FALSE),
    col = c(
      list(arm = arm_colours, timepoint = c(T1 = "grey85", T2 = "grey55", T3 = "grey25")),
      map(tracks, grey_ramp)
    ),
    na_col = "grey80", annotation_name_gp = gpar(fontsize = 7), annotation_name_side = "left",
    simple_anno_size = unit(3, "mm")
  )
  row_order <- paste(
    "Rows are ordered by average-linkage clustering (1 - Pearson r) of limpa's",
    "matrix, for display only."
  )
  effect_limit <- cap(contrast_results$effect)
  overviews <- imap(row_groups, \(rows, group) {
    m <- sample_values[rows, columns$sample_id, drop = FALSE]
    m <- sweep(m, 1, rowMeans(m, na.rm = TRUE))
    limit <- cap(m)
    labelled_heat <- \(...) heat(..., labels = feature_label)
    biopsy_map <- labelled_heat(m, "biopsies", c(-limit, 0, limit), "log2 abundance, centred",
      row_names = TRUE, column_title = character(0), top_annotation = sample_annotation,
      column_split = columns$group, cluster_column_slices = FALSE,
      column_labels = columns$subject, column_names_gp = gpar(fontsize = 5)
    )
    contrast_grid <- \(value) as_grid(contrast_results, rows, contrast_names, "contrast", value)
    contrast_maps <- labelled_heat(
      contrast_grid(\(d) d$effect), "log2 fold change",
      c(-effect_limit, 0, effect_limit), "log2 FC"
    ) +
      labelled_heat(contrast_grid(\(d) status_of(d$p, d$fdr)), "moderated t p and BH q",
        row_names = TRUE
      )
    n_rows <- format(length(rows), big.mark = ",")
    list(
      heatmap_page(biopsy_map, sprintf("Biopsy overview: %s", group), paste(
        sprintf("All %s %s by %d biopsies,", n_rows, feature_noun, ncol(m)),
        "grouped by arm and time point; column labels are subjects.",
        "Fill: log2 abundance minus the protein's mean over biopsies.",
        sprintf("Colours stop at the 99th percentile of |value|, %.2g.", limit),
        row_order, "Grey: missing."
      ), legends = list()),
      heatmap_page(contrast_maps, sprintf("Contrast overview: %s", group), paste(
        sprintf("All %s %s by the nine contrasts.", n_rows, feature_noun),
        "Fill: log2 fold change from the limma fit; status: p and BH q within contrast.",
        sprintf("Colours stop at the 99th percentile of |value|, %.2g.", effect_limit),
        row_order, "Grey: not tested. Table: c_data/02_contrasts.xlsx, DEP_matrix."
      ))
    )
  }) |>
    list_flatten()
  stopifnot(length(overviews) == 2 * length(row_groups))

  corrected <- readRDS(paths[["corrected"]])
  sample_meta <- as_tibble(meta)

  # The path of each arm's p < 0.05 training sets through T1, T2 and T3, in each arm: every
  # subject's set mean minus their own T1 value, then the mean and SE across subjects. A subject
  # with no T1 biopsy drops out.
  stopifnot(identical(rownames(corrected$data), rownames(abundance)))
  set_trajectory <- results |>
    filter(contrast %in% c("Training_C1", "Training_C2"), p < 0.05) |>
    mutate(set_from = contrast, set = if_else(logFC > 0, "nominal_up", "nominal_down")) |>
    summarise(members = list(uniprot_id), .by = c(set_from, set)) |>
    pmap(\(set_from, set, members) {
      tibble(
        sample_id = colnames(corrected$data),
        value = colMeans(corrected$data[members, , drop = FALSE])
      ) |>
        left_join(select(sample_meta, sample_id, subject, arm, timepoint), by = "sample_id") |>
        filter("T1" %in% timepoint, .by = subject) |>
        mutate(change = value - value[timepoint == "T1"], .by = subject) |>
        summarise(
          mean_change = mean(change), se = sd(change) / sqrt(n()), n_subjects = n(),
          .by = c(arm, timepoint)
        ) |>
        mutate(set_from = set_from, set = set, n_proteins = length(members), .before = 1)
    }) |>
    list_rbind()

  # Whole-proteome structure: PCA of the complete corrected matrix (proteins scaled), and PERMANOVA
  # on its Euclidean distances. Time is tested within subject (permutations stay inside each
  # person); arm is tested on one profile per subject (the mean of their samples), or on one time
  # point.
  scaled <- scale(t(corrected$data))
  pca <- prcomp(scaled)
  pca_variance <- tibble(
    component = paste0("PC", 1:4),
    variance_explained = round(summary(pca)$importance[2, 1:4], 3)
  )
  pca_scores <- as_tibble(pca$x[, 1:2], rownames = "sample_id") |>
    left_join(select(sample_meta, sample_id, subject, arm, timepoint), by = "sample_id")
  permanova <- function(rows, term, within_subject) {
    x <- scaled[rows, , drop = FALSE]
    data <- sample_meta[rows, ]
    if (!within_subject && anyDuplicated(data$subject)) {
      x <- rowsum(x, data$subject) / as.vector(table(data$subject))
      data <- distinct(data, subject, arm)[match(rownames(x), unique(data$subject)), ]
    }
    control <- if (within_subject) {
      permute::how(blocks = factor(data$subject), nperm = 9999)
    } else {
      permute::how(nperm = 9999)
    }
    set.seed(1)
    fit <- vegan::adonis2(dist(x) ~ ., data = data[term], permutations = control)
    tibble(n_samples = nrow(x), r2 = fit$R2[1], f = fit$F[1], p = fit$`Pr(>F)`[1])
  }
  time_pairs <- list(all = c("T1", "T2", "T3"), training = c("T1", "T2"), acute = c("T2", "T3"))
  global_tests <- bind_rows(
    imap(time_pairs, \(times, window) {
      map(c("both", "C1", "C2"), \(who) {
        rows <- sample_meta$timepoint %in% times & (who == "both" | sample_meta$arm == who)
        permanova(rows, "timepoint", within_subject = TRUE) |>
          mutate(term = "time", window = window, arms = who, .before = 1)
      })
    }) |>
      list_flatten() |>
      list_rbind(),
    map(c("subject mean", "T1", "T2", "T3"), \(when) {
      rows <- when == "subject mean" | sample_meta$timepoint == when
      permanova(rows, "arm", within_subject = FALSE) |>
        mutate(term = "arm", window = when, arms = "both", .before = 1)
    }) |>
      list_rbind()
  )
  # The same test on the phenotype: the six size and strength traits at T1 and T2, scaled, time
  # within subject. Arm is not tested, because C1 and C2 were split on how much these traits
  # changed.
  trait_sets <- phenotype_tables(trait_table)
  trait_rows <- filter(trait_sets$biopsy, if_all(all_of(trait_sets$traits), \(x) !is.na(x)))
  set.seed(1)
  trait_fit <- vegan::adonis2(
    dist(scale(trait_rows[trait_sets$traits])) ~ timepoint,
    data = trait_rows,
    permutations = permute::how(blocks = factor(trait_rows$subject), nperm = 9999)
  )
  global_tests <- global_tests |>
    mutate(data = "proteome", .before = 1) |>
    add_row(
      data = "phenotype", term = "time", window = "all", arms = "both",
      n_samples = nrow(trait_rows), r2 = trait_fit$R2[1], f = trait_fit$F[1],
      p = trait_fit$`Pr(>F)`[1]
    )
  pca_figure <- ggplot(pca_scores, aes(PC1, PC2, colour = timepoint, shape = arm)) +
    geom_point(size = 2) +
    scale_colour_manual(values = c(T1 = "#E69F00", T2 = "#0072B2", T3 = "#009E73")) +
    labs(
      title = "Samples on the complete, blood-corrected matrix",
      subtitle = with(
        filter(global_tests, data == "proteome", term == "time", window == "all", arms == "both"),
        sprintf("PERMANOVA, time within subject: R2 %.3f, p = %.4f", r2, p)
      ),
      x = sprintf("PC1 (%.1f%%)", 100 * pca_variance$variance_explained[1]),
      y = sprintf("PC2 (%.1f%%)", 100 * pca_variance$variance_explained[2]),
      caption = "Proteins scaled to unit variance. Table: c_data/02_contrasts.xlsx, pca, permanova."
    ) +
    theme_minimal(base_size = 9) +
    theme(plot.caption = element_text(hjust = 0))

  pdf(file.path(figures, "02_contrasts_figures.pdf"), width = 11, height = 8.5)
  walk(overviews, \(page) page$draw())
  walk(list(p_histogram, volcano, pca_figure), print)
  invisible(dev.off())

  write_workbook(
    list(
      DEP_matrix = results |>
        mutate(n_obs = n_obs[uniprot_id]) |>
        select(uniprot_id, gene, n_obs, contrast, logFC, ci_lower, ci_upper, p, fdr, qvalue) |>
        pivot_wider(
          names_from = contrast, values_from = c(logFC, ci_lower, ci_upper, p, fdr, qvalue),
          names_glue = "{contrast}_{.value}"
        ),
      contrast_summary = contrast_summary,
      null_calibration = null_calibration,
      detectable_effect = tibble(
        contrast = "Training_Interaction", quantile = c(0.5, 0.9),
        log2_fc = round(quantile(detectable, c(0.5, 0.9), na.rm = TRUE, names = FALSE), 3)
      ),
      fry_concordance = fry_concordance,
      set_trajectory = set_trajectory,
      unadjusted = unadjusted,
      pca = pca_scores,
      pca_variance = pca_variance,
      permanova = global_tests
    ),
    c(
      "logFC with its 95% CI, p, BH and Storey q per protein and contrast, with samples observed.",
      "Tested, p < 0.05 and BH counts per contrast, with its role.",
      "propTrueNull, median |t| and share of p below 0.2 per contrast.",
      "log2 effect the primary contrast detects at 80% power, median and 90th percentile.",
      paste(
        "fry: one arm's responders tested in the other arm's contrast, and training sets on each",
        "arm's acute contrast; with mixed p, mean log2 FC in both contrasts and sign kept."
      ),
      "Training p < 0.05 sets: mean change from each subject's T1 (and SE) per arm and time point.",
      "logFC and p per protein and contrast before blood adjustment, for sensitivity.",
      "PC1 and PC2 per sample, from the complete blood-corrected matrix with proteins scaled.",
      "Share of variance on the first four components.",
      paste(
        "PERMANOVA (9,999 permutations): proteome time per window and arm, arm per time point;",
        "phenotype time."
      )
    ),
    file.path(out, "02_contrasts.xlsx"),
    inputs, paths
  )
}
purrr::walk(inputs_used, run_input)
stack_inputs(here("03_Cluster_groups", "Proteome", "03_DEP"))

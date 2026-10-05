# The nine contrasts on limpa's quantification: limma with precision weights from limpa's standard
# errors, subject blocked, blood as a covariate; BH within each contrast, null checks, treat() at
# 1.15-fold, and fry on one arm's responders in the other arm.

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

stage <- here("02_Proteins", "02_Contrasts")
out <- file.path(stage, "c_data")
figures <- file.path(stage, "b_reports")
for (path in c(out, figures)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
inputs <- c(
  design = "02_Proteins/01_Design/c_data/design.rds",
  proteins = "01_Preprocess/03_Quantify/c_data/proteins.rds",
  corrected = "01_Preprocess/03_Quantify/c_data/DAList_corrected.rds",
  phenotype = "00_Input/phenotype.csv"
)
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
  topTable(eb, coef = contrast, number = Inf, sort.by = "none") |>
    as_tibble() |>
    select(logFC, AveExpr, t, P.Value, adj.P.Val) |>
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

ft <- treat(eb, lfc = log2(1.15), robust = TRUE)
treat_counts <- tibble(
  contrast = contrast_names,
  n_fdr_05 = colSums(decideTests(ft, p.value = 0.05) != 0, na.rm = TRUE)[contrast_names]
)

# fry on one arm's responders, tested in the other arm's contrast. Only mirrored within-arm pairs
# share no cell mean. fry reads limpa's complete matrix with the same design and block.
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
  "Training_HR", "Training_LR",
  "Training_LR", "Training_HR",
  "Acute_HR", "Acute_LR",
  "Acute_LR", "Acute_HR"
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
  fry(abundance,
    index = sets, design = design, contrast = contrast_matrix[, ranked_on],
    block = meta$subject, correlation = fry_correlation
  ) |>
    as_tibble(rownames = "set") |>
    transmute(
      set_from, ranked_on, set,
      n_proteins = NGenes, direction = Direction, p = PValue, fdr = FDR
    )
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
arm_colours <- c(HR = "#2166AC", LR = "#B2182B")
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
    rep(c("HR", "LR"), each = 3), c("T1", "T2", "T3"),
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

# ANCOVA, like the rank tests, reads limpa's complete blood-corrected matrix.
corrected <- readRDS(paths[["corrected"]])
ancova_results <- ancova(corrected$data, as_tibble(meta)) |>
  mutate(fdr = p.adjust(p, "BH"), qvalue = storey_q(p), .by = window)
# Leave one subject out on the tests that feed the candidate tables: the four HR-against-LR
# contrasts and both ANCOVA windows. A hit that stays at p < 0.05 in every rerun does not rest on
# one person.
sample_meta <- as_tibble(meta)
robustness <- leave_one_out(unique(sample_meta$subject), \(keep) {
  kept <- sample_meta$subject %in% keep
  bind_rows(
    fit_limpa(
      proteins[, kept], design[kept, ], contrast_matrix[, hr_lr_contrasts],
      sample_meta$subject[kept]
    ) |>
      contrast_rows(),
    ancova(observed_values[, kept], sample_meta[kept, ]) |>
      transmute(feature, test = paste("ANCOVA", window), p)
  ) |>
    select(feature, test, p)
})
ancova_figure <- ancova_results |>
  mutate(call = case_when(fdr < 0.05 ~ "BH < 0.05", p < 0.05 ~ "p < 0.05", TRUE ~ "not nominal")) |>
  ggplot(aes(effect, -log10(p), colour = call)) +
  geom_point(size = 0.8, alpha = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "grey50") +
  facet_wrap(~window) +
  scale_colour_manual(
    values = c(`BH < 0.05` = "#B2182B", `p < 0.05` = "#2166AC", `not nominal` = "grey75"),
    name = NULL
  ) +
  labs(
    title = "ANCOVA: HR minus LR at the later time point, at equal starting value",
    subtitle = sprintf(
      "Training: %d of %d proteins at p < 0.05. Acute: %d of %d.",
      sum(ancova_results$p[ancova_results$window == "training"] < 0.05),
      sum(ancova_results$window == "training"),
      sum(ancova_results$p[ancova_results$window == "acute"] < 0.05),
      sum(ancova_results$window == "acute")
    ),
    x = "log2 abundance, HR minus LR", y = "-log10 p",
    caption = paste(
      "lm(later ~ earlier + arm) per protein; training is T2 on T1, acute T3 on T2. Dashed line:",
      "p = 0.05. BH within window. Table: c_data/02_contrasts.xlsx, ancova."
    )
  ) +
  theme_minimal(base_size = 9) +
  theme(legend.position = "top", plot.caption = element_text(hjust = 0))

pdf(file.path(figures, "02_contrasts_figures.pdf"), width = 11, height = 8.5)
walk(overviews, \(page) page$draw())
walk(list(p_histogram, volcano, ancova_figure), print)
invisible(dev.off())

write_workbook(
  list(
    DEP_matrix = results |>
      mutate(n_obs = n_obs[uniprot_id]) |>
      select(uniprot_id, gene, n_obs, contrast, logFC, p, fdr, qvalue) |>
      pivot_wider(
        names_from = contrast, values_from = c(logFC, p, fdr, qvalue),
        names_glue = "{contrast}_{.value}"
      ),
    contrast_summary = contrast_summary,
    null_calibration = null_calibration,
    treat = treat_counts,
    detectable_effect = tibble(
      contrast = "Training_Interaction", quantile = c(0.5, 0.9),
      log2_fc = round(quantile(detectable, c(0.5, 0.9), na.rm = TRUE, names = FALSE), 3)
    ),
    fry_concordance = fry_concordance,
    unadjusted = unadjusted,
    ancova = ancova_results,
    leave_one_out = robustness
  ),
  c(
    "logFC, p, BH and Storey q per protein and contrast, with samples observed.",
    "Tested, p < 0.05 and BH counts per contrast, with its role.",
    "propTrueNull, median |t| and share of p below 0.2 per contrast.",
    "Proteins at BH < 0.05 under treat() with a 1.15-fold floor, per contrast.",
    "log2 effect the primary contrast detects at 80% power, median and 90th percentile.",
    "fry: one arm's responders tested in the other arm's contrast.",
    "logFC and p per protein and contrast before blood adjustment, for sensitivity.",
    "ANCOVA per protein and window: arm effect at equal starting value, t, p, BH FDR and Storey q.",
    "Leave one subject out: reruns, reruns at p < 0.05, largest p and its subject, per test."
  ),
  file.path(out, "02_contrasts.xlsx"),
  inputs, paths
)

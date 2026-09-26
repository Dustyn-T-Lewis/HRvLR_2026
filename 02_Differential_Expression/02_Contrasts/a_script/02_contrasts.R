# The nine contrasts on the unimputed protein matrix: proteoDA's limma fit, BH within each
# contrast, null checks, treat() at 1.15-fold, and fry on one arm's responders in the other arm.

suppressPackageStartupMessages({
  library(here)
  library(proteoDA)
  library(limma)
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(tibble)
  library(tidyr)
  library(writexl)
})

stage <- here("02_Differential_Expression", "02_Contrasts")
out <- file.path(stage, "c_data")
figures <- file.path(stage, "b_reports")
for (path in c(out, figures)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
inputs <- c(
  design = "02_Differential_Expression/01_Design/c_data/design.rds",
  imputed = "01_Preprocess/03_Imputation/c_data/DAList_imputed.rds"
)
paths <- map_chr(inputs, here)
stopifnot(file.exists(paths))
d <- readRDS(paths[["design"]])
dal <- d$dal
contrast_names <- colnames(dal$design$contrast_matrix)
stopifnot(identical(rownames(dal$design$design_matrix), colnames(dal$data)))

fit <- fit_limma_model(dal)
res <- extract_DA_results(fit, pval_thresh = 0.05, lfc_thresh = 0, adj_method = "BH")
eb <- fit$eBayes_fit
# write_limma_plots() knits its own R Markdown, so it runs in a fresh R session. Its interactive
# reports land in b_reports/ and git ignores them.
invisible(callr::r(\(res, figures) {
  proteoDA::write_limma_plots(res,
    grouping_column = "group", output_dir = figures,
    table_columns = c("uniprot_id", "gene", "protein"), title_column = "gene", overwrite = TRUE
  )
}, args = list(res = res, figures = figures)))

# Effect the primary contrast could detect at 80% power, two-sided 0.05.
se <- eb$stdev.unscaled[, "Training_Interaction"] * sqrt(eb$s2.post)
detectable <- se * (qt(0.975, eb$df.total) + qt(0.80, eb$df.total))

annotation <- fit$annotation[, c("uniprot_id", "gene", "protein", "description")]
stopifnot(identical(rownames(res$results[[1]]), rownames(annotation)))
# A protein with an empty cell in a contrast comes back NA there, so each contrast has its own
# denominator. pi_score (Xiao et al. 2014, inverted) only selects sets for the fry check below.
results <- res$results |>
  map(\(r) bind_cols(annotation, r)) |>
  list_rbind(names_to = "contrast") |>
  as_tibble() |>
  rename(p = P.Value, fdr = adj.P.Val) |>
  mutate(
    pi_score = p^abs(logFC),
    sig_pi = case_when(
      pi_score < 0.05 & logFC > 0 ~ 1L,
      pi_score < 0.05 & logFC < 0 ~ -1L,
      TRUE ~ 0L
    )
  )
n_obs <- rowSums(!is.na(fit$data))

contrast_summary <- results |>
  summarise(
    n_tested = sum(!is.na(p)),
    n_nominal = sum(p < 0.05, na.rm = TRUE),
    n_fdr_05 = sum(fdr < 0.05, na.rm = TRUE),
    n_fdr_10 = sum(fdr < 0.10, na.rm = TRUE),
    min_fdr = min(fdr, na.rm = TRUE),
    .by = contrast
  ) |>
  mutate(
    n_expected = 0.05 * n_tested, ratio = round(n_nominal / n_expected, 2), .after = n_nominal
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
# share no cell mean. fry takes no missing value, so it reads the imputed matrix with the same
# design and block.
imputed <- readRDS(paths[["imputed"]])
abundance <- as.matrix(imputed$data)
design <- fit$design$design_matrix
meta <- fit$metadata
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
    index = sets, design = design, contrast = fit$design$contrast_matrix[, ranked_on],
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
    title = "p-value distribution per contrast", x = "p", y = "proteins",
    caption = paste(
      "Dashed line: a flat null (tested / 20). A spike at zero on a flat background is signal;",
      "a slope up toward 1 is a conservative test. Table: c_data/02_contrasts.xlsx, DEP_matrix."
    )
  ) +
  theme_minimal(base_size = 10)

# Every protein nominal in at least one contrast, 75 rows to a page: fill is log2 FC, size is
# -log10 p, a black ring marks BH < 0.05.
label_of <- annotation |>
  mutate(label = if_else(
    is.na(gene) | duplicated(gene) | duplicated(gene, fromLast = TRUE),
    paste0(coalesce(gene, "?"), " (", uniprot_id, ")"), gene
  )) |>
  with(set_names(label, uniprot_id))
ranked <- results |>
  filter(!is.na(p)) |>
  summarise(n_nominal = sum(p < 0.05), best = min(p), .by = uniprot_id) |>
  filter(n_nominal > 0) |>
  arrange(desc(n_nominal), best) |>
  mutate(page = ceiling(row_number() / 75))
hit_pages <- map(split(ranked$uniprot_id, ranked$page), \(ids) {
  results |>
    filter(uniprot_id %in% ids) |>
    mutate(
      label = factor(label_of[uniprot_id], levels = rev(label_of[ids])),
      contrast = factor(contrast, levels = contrast_names),
      nominal = !is.na(p) & p < 0.05
    ) |>
    ggplot(aes(contrast, label)) +
    geom_point(data = \(x) filter(x, !nominal), colour = "grey85", size = 0.6) +
    geom_point(
      data = \(x) filter(x, nominal), aes(size = -log10(p), fill = logFC),
      shape = 21, colour = "grey30", stroke = 0.2
    ) +
    geom_point(
      data = \(x) filter(x, fdr < 0.05), aes(size = -log10(p)),
      shape = 21, colour = "black", stroke = 0.9
    ) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B") +
    scale_size_continuous(range = c(1, 4.5), limits = c(-log10(0.05), NA)) +
    scale_x_discrete(drop = FALSE) +
    labs(
      title = sprintf("Proteins nominal in at least one contrast (%d)", nrow(ranked)),
      x = NULL, y = NULL, fill = "log2 FC", size = "-log10 p",
      caption = paste(
        "Rows: proteins at p < 0.05 in at least one contrast, most contrasts first, then best p.",
        "Grey speck: tested, not nominal. Black ring: BH < 0.05. Chance: 5% of tested per",
        "contrast (contrast_summary). Table: c_data/02_contrasts.xlsx, DEP_matrix."
      )
    ) +
    theme_minimal(base_size = 9) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1), axis.text.y = element_text(size = 5.5)
    )
})

saveRDS(res, file.path(out, "fit.rds"), compress = "xz")
pdf(file.path(figures, "02_contrasts_figures.pdf"), width = 11, height = 8.5)
walk(c(list(p_histogram), hit_pages), print)
invisible(dev.off())

sheets <- list(
  DEP_matrix = results |>
    mutate(n_obs = n_obs[uniprot_id]) |>
    select(uniprot_id, gene, n_obs, contrast, logFC, p, fdr) |>
    pivot_wider(
      names_from = contrast, values_from = c(logFC, p, fdr), names_glue = "{contrast}_{.value}"
    ),
  contrast_summary = contrast_summary,
  null_calibration = null_calibration,
  treat = treat_counts,
  detectable_effect = tibble(
    contrast = "Training_Interaction", quantile = c(0.5, 0.9),
    log2_fc = round(quantile(detectable, c(0.5, 0.9), na.rm = TRUE, names = FALSE), 3)
  ),
  fry_concordance = fry_concordance,
  input_manifest = tibble(
    input = names(inputs), path = unname(inputs), md5 = unname(tools::md5sum(paths))
  ),
  package_versions = sessioninfo::package_info("loaded", dependencies = FALSE) |>
    as_tibble() |>
    select(package, version = loadedversion, source)
)
read_me <- tibble(sheet = names(sheets), holds = c(
  "logFC, p and BH FDR per protein and contrast, with samples observed.",
  "Tested, nominal, chance-expected and BH counts per contrast, with its role.",
  "propTrueNull, median |t| and share of p below 0.2 per contrast.",
  "Proteins at BH < 0.05 under treat() with a 1.15-fold floor, per contrast.",
  "log2 effect the primary contrast detects at 80% power, median and 90th percentile.",
  "fry: one arm's responders tested in the other arm's contrast.",
  "Files read, with md5.",
  "Packages loaded at run time."
))
write_xlsx(c(list(read_me = read_me), sheets), file.path(out, "02_contrasts.xlsx"))

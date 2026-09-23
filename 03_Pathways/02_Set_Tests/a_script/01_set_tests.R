# Set-level tests of the nine contrasts, two ways.
#
# fry is the inferential test. It rotates residuals under the same design,
# subject block and within-subject correlation as the protein fit, so it does
# not assume proteins are independent; METHOD_RANKING.md ranks it first for
# that reason. fry cannot take an NA, so it runs on the missForest matrix; the
# design and contrasts are the fit's own.
#
# fgsea on the moderated t is the display layer. It supplies the NES and
# leading edge the packet draws, and collapsePathways marks which significant
# sets are not just restatements of a larger one. Its gene-permutation null
# treats proteins as exchangeable, so its p is not the one to quote.

pacman::p_load(here, dplyr, tidyr, tibble, purrr, limma, fgsea)

source(here("functions", "contrasts.R"))
source(here("functions", "shared_pathway_utils.R"))
source(here("functions", "shared_utils.R"))

set.seed(42)

OUT_DIR <- here("03_Pathways", "02_Set_Tests", "c_data")

gs <- readRDS(here("03_Pathways", "01_Gene_Sets", "c_data", "gene_sets.rds"))
proteins <- readRDS(here(
  "02_Proteins", "01_Differential", "c_data", "proteins.rds"
))
fit <- readRDS(here(
  "02_Proteins", "01_Differential", "c_data", "01_limma_DAList.rds"
))
imputed <- readRDS(here(
  "01_Preprocess", "03_Imputation", "c_data", "DAList_imputed_missforest.rds"
))

abund <- as.matrix(imputed$data)
design <- fit$design$design_matrix
contrasts <- fit$design$contrast_matrix
subject <- fit$metadata$subject
stopifnot(
  identical(colnames(abund), rownames(design)),
  identical(rownames(abund), rownames(proteins$abund)),
  !anyNA(abund)
)
correlation <- duplicateCorrelation(abund, design, block = subject)$consensus
index <- map(gs$sets, \(ids) match(ids, rownames(abund)))

fry_res <- map(CONTRAST_NAMES, function(ct) {
  fry(abund,
    index = index, design = design, contrast = contrasts[, ct],
    block = subject, correlation = correlation, sort = "none"
  ) |>
    as_tibble(rownames = "set") |>
    transmute(
      contrast = ct, set = .data$set, direction = .data$Direction,
      fry_p = .data$PValue
    )
}) |>
  list_rbind()

fgsea_res <- map(CONTRAST_NAMES, function(ct) {
  stats <- proteins$stats$t[, ct]
  stats <- sort(stats[!is.na(stats)], decreasing = TRUE)
  res <- fgsea(gs$sets,
    stats = stats, minSize = SET_FLOOR, maxSize = SET_CEILING
  )
  hits <- res[res$padj < 0.05, ]
  main <- if (nrow(hits)) {
    collapsePathways(hits, gs$sets, stats, pval.threshold = 0.05)$mainPathways
  } else {
    character()
  }
  as_tibble(res) |>
    transmute(
      contrast = ct, set = .data$pathway, nes = .data$NES,
      fgsea_p = .data$pval, fgsea_padj = .data$padj,
      leading_edge = map_chr(.data$leadingEdge, paste, collapse = ";"),
      main = .data$pathway %in% main
    )
}) |>
  list_rbind()

# BH within each contrast and collection, never across: 33 Hallmark sets and
# 1156 GO:BP sets are different-sized questions with different nulls.
set_tests <- gs$catalog |>
  select("collection", "set", "theme", "size_detected") |>
  inner_join(fry_res, by = "set", relationship = "one-to-many") |>
  left_join(fgsea_res, by = c("contrast", "set")) |>
  mutate(
    contrast = factor(.data$contrast, levels = CONTRAST_NAMES),
    fry_fdr = p.adjust(.data$fry_p, "BH"),
    .by = c("contrast", "collection")
  )

set_matrix <- function(col) {
  set_tests |>
    select("set", "contrast", value = all_of(col)) |>
    pivot_wider(names_from = "contrast", values_from = "value") |>
    column_to_rownames("set") |>
    as.matrix()
}
signed_fry <- set_tests |>
  mutate(value = -log10(.data$fry_p) * if_else(.data$direction == "Up", 1, -1))

clear_dir(OUT_DIR)
saveRDS(
  list(
    set_tests = set_tests,
    nes = set_matrix("nes"),
    fry_fdr = set_matrix("fry_fdr"),
    fry_signed = signed_fry |>
      select("set", "contrast", "value") |>
      pivot_wider(names_from = "contrast", values_from = "value") |>
      column_to_rownames("set") |>
      as.matrix(),
    correlation = correlation
  ),
  file.path(OUT_DIR, "set_tests.rds")
)
openxlsx::write.xlsx(c(
  list(summary = set_tests |>
    summarise(
      n_sets = n(),
      fry_nominal = sum(.data$fry_p < 0.05),
      fry_expected = 0.05 * n(),
      fry_fdr05 = sum(.data$fry_fdr < 0.05),
      fgsea_padj05 = sum(.data$fgsea_padj < 0.05, na.rm = TRUE),
      fgsea_main = sum(.data$main, na.rm = TRUE),
      .by = c("contrast", "collection")
    )),
  split(set_tests, set_tests$contrast)
), file.path(OUT_DIR, "set_tests.xlsx"))

message(sprintf("within-subject correlation, imputed: %.3f", correlation))
print(
  count(
    filter(set_tests, .data$fry_fdr < 0.05), .data$contrast, .data$collection
  ),
  n = Inf
)

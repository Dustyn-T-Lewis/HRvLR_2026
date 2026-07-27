#!/usr/bin/env Rscript
# The proteins the model never tests, and the true BH denominator.
#
# min_groups = 1 admits a protein detected in a single Group_Time cell, on the reasoning that a
# protein switched on in one responder group and off in the other is the biology under study.
# With ~ 0 + group cell-means coding, a cell with no observations has no estimable mean, so limma
# returns NA for every contrast touching it. Those proteins are admitted by the filter and then
# silently never tested: the very on/off signal the filter was widened to keep is invisible to
# the model that follows.
#
# Nothing here changes the fit. It reports what the fit could not reach, so the gap between the
# filter's intent and the model's reach is on the record rather than hidden in NA rows.

pacman::p_load(proteoDA, here, readr, dplyr, tibble, purrr, tidyr)

dal <- readRDS(here("01_Filtering", "c_data", "DAList_filtered.rds"))
fit <- readRDS(here("03_DEP", "a_non_imputed", "c_data", "01_limma_DAList.rds"))

cells <- factor(dal$metadata$Group_Time[match(colnames(dal$data), dal$metadata$Col_ID)])
observed <- t(apply(dal$data, 1, \(x) tapply(!is.na(x), cells, sum)))

untested <- tibble(
  uniprot_id = rownames(observed),
  n_empty_cells = rowSums(observed == 0),
  n_detected = rowSums(!is.na(dal$data))
) |>
  filter(n_empty_cells > 0) |>
  left_join(as_tibble(fit$annotation), by = "uniprot_id") |>
  bind_cols(as_tibble(observed[rowSums(observed == 0) > 0, , drop = FALSE])) |>
  select(uniprot_id, gene, protein, n_detected, n_empty_cells, HR_T1:LR_T3) |>
  arrange(desc(n_empty_cells), n_detected)

# The BH denominator is the non-NA count, not every row in the matrix.
denominators <- imap_dfr(fit$results, \(r, cname) {
  tibble(
    contrast = cname,
    n_tested = sum(!is.na(r$adj.P.Val)),
    n_untested = sum(is.na(r$adj.P.Val))
  )
})

out_dir <- here("03_DEP", "a_non_imputed", "b_reports")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write_csv(untested, file.path(out_dir, "untested_proteins.csv"))
write_csv(denominators, file.path(out_dir, "bh_denominators.csv"))

cat(sprintf(
  "%d proteins reach the model with >=1 empty Group_Time cell and are never tested.\n",
  nrow(untested)
))
print(count(untested, n_empty_cells, name = "n_proteins"), n = Inf)
cat(sprintf(
  "\nTested-N per contrast (the BH denominator; matrix has %d rows):\n",
  nrow(fit$data)
))
print(as.data.frame(denominators), row.names = FALSE)

#!/usr/bin/env Rscript
# The nine contrasts refitted through proteoDA with blood in the design.
#
# 05_blood_adjusted.R adjusts downstream, inside fit_feature_contrasts(). This
# one puts the covariate in the proteoDA model itself, so the moderated variance,
# the duplicateCorrelation consensus and the BH denominators are all estimated
# against the adjusted residual rather than borrowed from the primary fit.
#
# One covariate, not a panel of blood proteins. At 45 samples the design already
# spends six degrees of freedom on the group means; a second blood term buys
# almost nothing and costs power everywhere. The index is the mean log2 of five
# haemoglobins, and those five were removed by the contaminant filter, so it is
# built from proteins that are not among the 1900 being tested. Regressing the
# matrix on it is therefore not circular.
#
# The primary fit is untouched and stays the headline. This writes to its own
# directory and reports counts beside the primary, never promoting a row on its
# own. A protein significant only in a secondary model, in a study whose
# interactions are null, is a candidate and not a result.

pacman::p_load(proteoDA, here, dplyr, tibble, readr, openxlsx)
source(here("03_Features", "contrasts.R"))
source(here("functions", "blood_index_model.R"))
source(here("functions", "shared_utils.R"))

OUT <- here("03_Features", "01_Proteins", "c_data", "blood_covariate")
clear_dir(OUT)

dal <- readRDS(here("02_Normalization", "c_data", "DAList_normalized.rds"))
blood <- blood_index_data()

meta <- as.data.frame(dal$metadata)
meta$group <- factor(meta$Group_Time, levels = GROUP_LEVELS)
meta$subject <- meta$Subject_ID
meta$blood <- blood$blood_index[match(meta$Col_ID, blood$Col_ID)]
stopifnot(!anyNA(meta$blood))
dal$metadata <- meta

# add_design() strips the variable name off each design column, which is what
# turns groupHR_T1 into HR_T1. A continuous term has nothing left after that
# strip, so the covariate arrives unnamed and makeContrasts refuses the matrix.
# Naming it back is the whole fix; the column itself is correct.
name_covariate <- function(d, nm) {
  cols <- colnames(d$design$design_matrix)
  blank <- !nzchar(cols)
  if (any(blank)) {
    cols[blank] <- nm
    colnames(d$design$design_matrix) <- cols
  }
  d
}

fit_with <- function(formula) {
  d <- name_covariate(add_design(dal, formula), "blood")
  d <- add_contrasts(d, contrasts_vector = HRVLR_CONTRASTS)
  extract_DA_results(
    fit_limma_model(d),
    pval_thresh = 0.10, lfc_thresh = 0, adj_method = "BH"
  )
}

primary <- fit_with("~ 0 + group + (1 | subject)")
adjusted <- fit_with("~ 0 + group + blood + (1 | subject)")

as_long <- function(res, tag) {
  bind_rows(lapply(names(res$results), function(ct) {
    as_tibble(res$results[[ct]], rownames = "uniprot_id") |>
      transmute(
        contrast = ct, uniprot_id,
        logFC = .data$logFC, p = .data$P.Value, bh = .data$adj.P.Val
      )
  })) |>
    rename_with(\(x) paste0(x, "_", tag), c("logFC", "p", "bh"))
}

ann <- as.data.frame(dal$annotation)
compare <- as_long(primary, "primary") |>
  inner_join(as_long(adjusted, "adj"), by = c("contrast", "uniprot_id")) |>
  mutate(gene = ann$gene[match(uniprot_id, ann$uniprot_id)])

summary_tbl <- compare |>
  summarise(
    n = dplyr::n(),
    bh05_primary = sum(bh_primary < 0.05, na.rm = TRUE),
    bh05_adj = sum(bh_adj < 0.05, na.rm = TRUE),
    bh10_adj = sum(bh_adj < 0.10, na.rm = TRUE),
    nominal_primary = sum(p_primary < 0.05, na.rm = TRUE),
    nominal_adj = sum(p_adj < 0.05, na.rm = TRUE),
    min_bh_adj = min(bh_adj, na.rm = TRUE),
    .by = contrast
  ) |>
  mutate(responder = contrast %in% RESPONDER_CONTRASTS) |>
  arrange(min_bh_adj)

print(as.data.frame(summary_tbl), row.names = FALSE, digits = 3)

gained <- compare |>
  filter(bh_adj < 0.05, bh_primary >= 0.05) |>
  arrange(bh_adj)
cat(sprintf(
  "\nproteins reaching BH < .05 only after adjustment: %d (%d in a responder contrast)\n",
  nrow(gained), sum(gained$contrast %in% RESPONDER_CONTRASTS)
))
if (nrow(gained)) {
  gained |>
    transmute(contrast, gene,
      lfc = round(logFC_adj, 2),
      p = signif(p_adj, 2), bh = round(bh_adj, 4)
    ) |>
    head(20) |>
    as.data.frame() |>
    print(row.names = FALSE)
}

write.xlsx(
  list(summary = summary_tbl, gained = gained, all = compare),
  file.path(OUT, "blood_covariate_contrasts.xlsx")
)
cat(sprintf("\nwrote %s\n", file.path(OUT, "blood_covariate_contrasts.xlsx")))

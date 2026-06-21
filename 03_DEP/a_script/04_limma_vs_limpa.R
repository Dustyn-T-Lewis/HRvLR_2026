# HRvLR — Head-to-head: limma+missForest vs LIMPA (DPC + voomaLmFit)
# Per contrast: DEP counts, set overlap, logFC concordance, top disagreements.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(purrr)
})

setwd(rprojroot::find_rstudio_root_file())

cfg <- list(
  limma_dir = "03_DEP/c_data/04_per_contrast_results",
  limpa_dir = "03_DEP/c_data/limpa/04_per_contrast_results",
  out_dir   = "03_DEP/c_data/limpa_compare",
  contrasts = c("Training_HR", "Training_LR", "Acute_HR", "Acute_LR",
                "Baseline_HRvLR", "Trained_HRvLR", "Acute_HRvLR",
                "Training_Interaction", "Acute_Interaction")
)
dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)

read_one <- function(dir, cn) {
  read_csv(file.path(dir, paste0(cn, ".csv")), show_col_types = FALSE) |>
    select(uniprot_id, gene, logFC, P.Value, adj.P.Val, pi_score)
}

compare_contrast <- function(cn) {
  lm <- read_one(cfg$limma_dir, cn) |>
    rename(logFC_lm = logFC, P_lm = P.Value, FDR_lm = adj.P.Val, pi_lm = pi_score)
  lp <- read_one(cfg$limpa_dir, cn) |>
    rename(logFC_lp = logFC, P_lp = P.Value, FDR_lp = adj.P.Val, pi_lp = pi_score) |>
    select(-gene)

  j <- inner_join(lm, lp, by = "uniprot_id")
  s <- function(x) sum(x, na.rm = TRUE)

  summary_row <- tibble(
    contrast            = cn,
    n_proteins          = nrow(j),
    fdr05_lm            = s(j$FDR_lm < 0.05),
    fdr05_lp            = s(j$FDR_lp < 0.05),
    fdr05_overlap       = s(j$FDR_lm < 0.05 & j$FDR_lp < 0.05),
    fdr05_jaccard       = round(fdr05_overlap /
      max(1, s(j$FDR_lm < 0.05 | j$FDR_lp < 0.05)), 3),
    fdr10_lm            = s(j$FDR_lm < 0.10),
    fdr10_lp            = s(j$FDR_lp < 0.10),
    fdr10_overlap       = s(j$FDR_lm < 0.10 & j$FDR_lp < 0.10),
    fdr10_jaccard       = round(fdr10_overlap /
      max(1, s(j$FDR_lm < 0.10 | j$FDR_lp < 0.10)), 3),
    pi05_lm             = s(j$pi_lm < 0.05),
    pi05_lp             = s(j$pi_lp < 0.05),
    pi05_overlap        = s(j$pi_lm < 0.05 & j$pi_lp < 0.05),
    pi05_jaccard        = round(pi05_overlap /
      max(1, s(j$pi_lm < 0.05 | j$pi_lp < 0.05)), 3),
    spearman_logFC      = round(cor(j$logFC_lm, j$logFC_lp, method = "spearman",
                                    use = "complete.obs"), 3),
    pearson_logFC       = round(cor(j$logFC_lm, j$logFC_lp, method = "pearson",
                                    use = "complete.obs"), 3),
    median_abs_dlogFC   = round(median(abs(j$logFC_lm - j$logFC_lp),
                                       na.rm = TRUE), 3),
    sign_concordance    = round(mean(sign(j$logFC_lm) == sign(j$logFC_lp),
                                     na.rm = TRUE), 3)
  )

  disagree <- j |>
    mutate(
      flag = case_when(
        FDR_lm < 0.10 & FDR_lp >= 0.10 ~ "limma_only",
        FDR_lp < 0.10 & FDR_lm >= 0.10 ~ "limpa_only",
        TRUE ~ NA_character_
      )
    ) |>
    filter(!is.na(flag)) |>
    mutate(score = pmin(FDR_lm, FDR_lp, na.rm = TRUE)) |>
    arrange(score) |>
    head(20) |>
    transmute(contrast = cn, flag, uniprot_id, gene,
              logFC_lm, FDR_lm, pi_lm,
              logFC_lp, FDR_lp, pi_lp)

  list(summary = summary_row, detail = j |> mutate(contrast = cn),
       disagree = disagree)
}

all_contrasts <- lapply(cfg$contrasts, compare_contrast)
summary_tbl   <- map_dfr(all_contrasts, "summary")
disagree_tbl  <- map_dfr(all_contrasts, "disagree")
detail_tbl    <- map_dfr(all_contrasts, "detail")

write_csv(summary_tbl,  file.path(cfg$out_dir, "01_compare_summary.csv"))
write_csv(disagree_tbl, file.path(cfg$out_dir, "02_top_disagreements.csv"))
write_csv(detail_tbl,   file.path(cfg$out_dir, "03_per_protein_detail.csv"))

cat("\n=== HRvLR limma vs LIMPA — summary ===\n")
print(summary_tbl |>
        select(contrast, n_proteins,
               fdr10_lm, fdr10_lp, fdr10_overlap, fdr10_jaccard,
               pi05_lm, pi05_lp, pi05_overlap, pi05_jaccard,
               spearman_logFC, sign_concordance) |>
        as.data.frame())

cat(sprintf("\nWritten to %s\n", cfg$out_dir))

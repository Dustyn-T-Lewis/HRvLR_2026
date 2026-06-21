# HRvLR DEP — LIMPA alternative (Path B: protein-level DPC + voomaLmFit)
# Runs alongside 01_run_dep.R (limma+missForest) for head-to-head comparison.
# Outputs to c_data/limpa/ to leave limma artifacts intact.
#
# 9 contrasts (Training/Acute within HR & LR, baseline/trained/acute HRvLR,
# Training & Acute interactions). See 01_run_dep.R for definitions.

suppressPackageStartupMessages({
  library(limpa)
  library(limma)
  library(dplyr)
  library(readr)
  library(tibble)
})

setwd(rprojroot::find_rstudio_root_file())
set.seed(42)

cfg <- list(
  norm_csv    = "01_normalization/c_data/02_normalized.csv",
  norm_rds    = "01_normalization/c_data/03_DAList_normalized.rds",
  out_dir     = "03_DEP/c_data/limpa",
  per_dir     = "03_DEP/c_data/limpa/04_per_contrast_results",
  pval_thresh = 0.10,
  pi_thresh   = 0.05,
  adj_method  = "BH",
  levels      = c("HR_T1", "HR_T2", "HR_T3", "LR_T1", "LR_T2", "LR_T3"),
  contrasts   = c(
    Training_HR          = "HR_T2 - HR_T1",
    Training_LR          = "LR_T2 - LR_T1",
    Acute_HR             = "HR_T3 - HR_T2",
    Acute_LR             = "LR_T3 - LR_T2",
    Baseline_HRvLR       = "HR_T1 - LR_T1",
    Trained_HRvLR        = "HR_T2 - LR_T2",
    Acute_HRvLR          = "HR_T3 - LR_T3",
    Training_Interaction = "(HR_T2 - HR_T1) - (LR_T2 - LR_T1)",
    Acute_Interaction    = "(HR_T3 - HR_T2) - (LR_T3 - LR_T2)"
  )
)

dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$per_dir, recursive = TRUE, showWarnings = FALSE)

# --- Load inputs ---
df <- read_csv(cfg$norm_csv, show_col_types = FALSE)
ann_cols <- c("uniprot_id", "protein", "gene", "description")
ann      <- df[, ann_cols]
mat      <- as.matrix(df[, setdiff(names(df), ann_cols)])
rownames(mat) <- ann$uniprot_id

dal_norm <- readRDS(cfg$norm_rds)
meta     <- as.data.frame(dal_norm$metadata)
mat      <- mat[, meta$Col_ID]

cat(sprintf("Loaded: %d proteins x %d samples | NA=%.1f%%\n",
            nrow(mat), ncol(mat), 100 * mean(is.na(mat))))
print(table(meta$Group, meta$Timepoint))

# --- LIMPA fit ---
y <- new("EList", list(E = mat,
                       genes = data.frame(uniprot_id = rownames(mat)),
                       targets = meta))

dpcfit <- dpc(y$E, verbose = FALSE)
cat(sprintf("DPC: intercept=%.3f slope=%.3f\n", dpcfit$dpc[1], dpcfit$dpc[2]))

imp_mask <- is.na(y$E)
y_imp    <- dpcQuantByRow(y, dpc = dpcfit, verbose = FALSE)
cat(sprintf("Imputed %d cells (%.1f%%)\n", sum(imp_mask), 100 * mean(imp_mask)))

group  <- factor(meta$Group_Time, levels = cfg$levels)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

fit <- voomaLmFitWithImputation(
  y_imp, design = design, imputed = imp_mask,
  block = meta$Subject_ID, sample.weights = TRUE, plot = FALSE
)
cm   <- do.call(makeContrasts, c(as.list(cfg$contrasts), list(levels = design)))
fit2 <- eBayes(contrasts.fit(fit, cm))

# --- Extract per-contrast results ---
results <- lapply(colnames(cm), function(cn) {
  tt <- topTable(fit2, coef = cn, number = Inf,
                 adjust.method = cfg$adj_method, sort.by = "none")
  tibble(uniprot_id = rownames(tt),
         logFC      = tt$logFC,
         P.Value    = tt$P.Value,
         adj.P.Val  = tt$adj.P.Val,
         t          = tt$t,
         B          = tt$B) |>
    mutate(pi_score = P.Value ^ abs(logFC),
           sig_pi   = case_when(
             pi_score < cfg$pi_thresh & logFC > 0 ~  1L,
             pi_score < cfg$pi_thresh & logFC < 0 ~ -1L,
             TRUE ~ 0L),
           contrast = cn)
}) |> setNames(colnames(cm))

# --- Per-contrast CSVs ---
for (cn in names(results)) {
  res <- results[[cn]] |>
    left_join(ann, by = "uniprot_id") |>
    select(uniprot_id, protein, gene, description, everything())
  write_csv(res, file.path(cfg$per_dir, paste0(cn, ".csv")))
}

# --- Combined wide table ---
combined <- ann
for (cn in names(results)) {
  r <- results[[cn]] |>
    select(uniprot_id, logFC, P.Value, adj.P.Val, pi_score, sig_pi) |>
    rename_with(~ paste0(.x, "_", cn), -uniprot_id)
  combined <- left_join(combined, r, by = "uniprot_id")
}
write_csv(combined, file.path(cfg$out_dir, "03_combined_results.csv"))

# --- Summary ---
da_summary <- bind_rows(lapply(names(results), function(cn) {
  r  <- results[[cn]]
  up <- r$logFC > 0; dn <- r$logFC < 0
  s  <- function(x) sum(x, na.rm = TRUE)
  tibble(
    contrast = cn,
    type     = c("up", "down", "nonsig"),
    sig.PVal = c(s(r$P.Value < 0.10 & up), s(r$P.Value < 0.10 & dn), s(r$P.Value >= 0.10)),
    sig.FDR  = c(s(r$adj.P.Val < 0.10 & up), s(r$adj.P.Val < 0.10 & dn), s(r$adj.P.Val >= 0.10)),
    sig.Pi   = c(s(r$sig_pi == 1L), s(r$sig_pi == -1L), s(r$sig_pi == 0L)),
    sig.FDR.05 = c(s(r$adj.P.Val < 0.05 & up), s(r$adj.P.Val < 0.05 & dn), s(r$adj.P.Val >= 0.05)),
    sig.FDR.10 = c(s(r$adj.P.Val < 0.10 & up), s(r$adj.P.Val < 0.10 & dn), s(r$adj.P.Val >= 0.10))
  )
}))
write_csv(da_summary, file.path(cfg$out_dir, "02_DA_summary.csv"))

# --- Save lightweight fit ---
saveRDS(list(results  = results,
             coef     = fit2$coefficients,
             dpc      = dpcfit$dpc,
             imp_frac = mean(imp_mask)),
        file.path(cfg$out_dir, "01_limpa_fit.rds"))

print(da_summary)
cat(sprintf("Done -> %s\n", cfg$out_dir))

writeLines(capture.output(sessionInfo()),
           file.path(cfg$out_dir, "sessionInfo.txt"))

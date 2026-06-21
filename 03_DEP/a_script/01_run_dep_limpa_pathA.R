# HRvLR DEP — LIMPA Path A: peptide-level pipeline (full LIMPA)
# Status: SCAFFOLDING — runs only when DIA-NN peptide-level output is present.
#
# When report from Auburn arrives:
#   1. Place DIA-NN report at: 00_input/dia_nn/report.tsv  (or report.parquet)
#   2. Create run-to-sample mapping at: 00_input/dia_nn/run_to_sample.csv
#      with columns: Run, Col_ID  (one row per LC-MS run; tech reps included)
#   3. Run this script. Outputs to c_data/limpa_pathA/.

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
  diann_report  = "00_input/dia_nn/report.tsv",
  diann_parquet = "00_input/dia_nn/report.parquet",
  run_map_csv   = "00_input/dia_nn/run_to_sample.csv",
  norm_rds      = "01_normalization/c_data/03_DAList_normalized.rds",
  out_dir       = "03_DEP/c_data/limpa_pathA",
  per_dir       = "03_DEP/c_data/limpa_pathA/04_per_contrast_results",
  q_cutoff      = 0.01,
  pi_thresh     = 0.05,
  pval_thresh   = 0.10,
  adj_method    = "BH",
  levels        = c("HR_T1", "HR_T2", "HR_T3", "LR_T1", "LR_T2", "LR_T3"),
  contrasts     = c(
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

# --- Pre-flight ---
report_path <- if (file.exists(cfg$diann_parquet)) cfg$diann_parquet else cfg$diann_report
if (!file.exists(report_path)) {
  stop(sprintf(paste0(
    "DIA-NN peptide-level report not found.\n",
    "  Expected at: %s (or %s)\n",
    "  Action: request `report.tsv` (or .parquet) from Auburn MS Core ",
    "and place at the path above.\n",
    "  Until then, use 01_run_dep_limpa.R (Path B, protein-level)."),
    cfg$diann_report, cfg$diann_parquet))
}
if (!file.exists(cfg$run_map_csv)) {
  stop(sprintf(paste0(
    "Run-to-sample mapping not found at %s.\n",
    "  Required columns: Run, Col_ID."), cfg$run_map_csv))
}

dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$per_dir, recursive = TRUE, showWarnings = FALSE)

# --- Read DIA-NN ---
y <- readDIANN(file = report_path, q.cutoffs = cfg$q_cutoff,
               log = TRUE, verbose = TRUE)
cat(sprintf("readDIANN: %d precursors x %d runs\n", nrow(y), ncol(y)))

run_map   <- read_csv(cfg$run_map_csv, show_col_types = FALSE)
dal_norm  <- readRDS(cfg$norm_rds)
norm_meta <- as.data.frame(dal_norm$metadata)

stopifnot(all(c("Run", "Col_ID") %in% colnames(run_map)))
missing_runs <- setdiff(colnames(y), run_map$Run)
if (length(missing_runs)) {
  stop("Runs in report not in run_map: ", paste(head(missing_runs, 5), collapse = ", "))
}

run_targets <- tibble(Run = colnames(y)) |>
  left_join(run_map, by = "Run") |>
  left_join(norm_meta, by = "Col_ID")
y$targets <- as.data.frame(run_targets)

keep <- !is.na(y$targets$Group_Time)
y <- y[, keep]
cat(sprintf("After QC: %d precursors x %d runs (%d samples)\n",
            nrow(y), ncol(y), length(unique(y$targets$Col_ID))))

# --- Peptide filters ---
y <- filterNonProteotypicPeptides(y)
y <- filterByDetection(y, n.samples = 3, n.detections = 1)
y <- filterSingletonPeptides(y)
y <- filterCompoundProteins(y)
cat(sprintf("After peptide filters: %d precursors\n", nrow(y)))

# --- DPC + protein-level summarization with SEs ---
dpcfit <- dpc(y$E, verbose = FALSE)
cat(sprintf("Peptide-level DPC: intercept=%.3f slope=%.3f\n",
            dpcfit$dpc[1], dpcfit$dpc[2]))

y_prot <- dpcQuant(y, protein.id = "Protein.Group", dpc = dpcfit, verbose = FALSE)
cat(sprintf("dpcQuant: %d proteins x %d runs\n", nrow(y_prot), ncol(y_prot)))

# --- Fit ---
meta   <- y_prot$targets
group  <- factor(meta$Group_Time, levels = cfg$levels)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

fit <- voomaLmFitWithImputation(
  y_prot, design = design, block = meta$Subject_ID,
  sample.weights = TRUE, plot = FALSE
)
cm   <- do.call(makeContrasts, c(as.list(cfg$contrasts), list(levels = design)))
fit2 <- eBayes(contrasts.fit(fit, cm))

results <- lapply(colnames(cm), function(cn) {
  tt <- topTable(fit2, coef = cn, number = Inf,
                 adjust.method = cfg$adj_method, sort.by = "none")
  tibble(uniprot_id = rownames(tt),
         logFC = tt$logFC, P.Value = tt$P.Value,
         adj.P.Val = tt$adj.P.Val, t = tt$t, B = tt$B) |>
    mutate(pi_score = P.Value ^ abs(logFC),
           sig_pi = case_when(
             pi_score < cfg$pi_thresh & logFC > 0 ~  1L,
             pi_score < cfg$pi_thresh & logFC < 0 ~ -1L,
             TRUE ~ 0L),
           contrast = cn)
}) |> setNames(colnames(cm))

for (cn in names(results)) {
  write_csv(results[[cn]], file.path(cfg$per_dir, paste0(cn, ".csv")))
}
saveRDS(list(results = results, coef = fit2$coefficients, dpc = dpcfit$dpc),
        file.path(cfg$out_dir, "01_limpa_pathA_fit.rds"))

cat(sprintf("Done -> %s\n", cfg$out_dir))
writeLines(capture.output(sessionInfo()),
           file.path(cfg$out_dir, "sessionInfo.txt"))

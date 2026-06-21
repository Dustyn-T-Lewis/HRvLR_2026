# compare/02_downstream.R — FC correlation, NES coherence, DEP counts
# Reads imp_list from parent environment (or CACHE_RDS)
# Writes: 02_Imputation/c_data/benchmark/02_downstream.csv

if (!exists("imp_list")) {
  source("02_Imputation/a_script/benchmark/_common.R")
  imp_list <- readRDS(CACHE_RDS)
}

suppressPackageStartupMessages({
  library(limma)
  library(proteoDA)
  library(fgsea)
  library(msigdbr)
})
select <- dplyr::select

# --- Build gene sets for NES comparison ---
# GO Slim + Hallmark (same as main pipeline)
hallmark <- msigdbr(species = "Homo sapiens", collection = "H") |>
  select(gs_name, gene_symbol) |>
  split(~gs_name) |>
  lapply(function(x) x$gene_symbol)

go_bp <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP") |>
  select(gs_name, gene_symbol) |>
  split(~gs_name) |>
  lapply(function(x) x$gene_symbol)

# Filter to moderate-size sets
gene_sets <- c(hallmark, go_bp)
gene_sets <- gene_sets[sapply(gene_sets, length) >= 15 & sapply(gene_sets, length) <= 500]
cat(sprintf("Using %d gene sets for NES comparison\n", length(gene_sets)))

# --- Helper: run limma on a matrix ---
# HRvLR model: Acute_LR contrast (strongest signal)
run_limma_benchmark <- function(mat, meta_df) {
  meta_df$group <- factor(meta_df$Group_Time,
    levels = c("HR_T1", "HR_T2", "HR_T3", "LR_T1", "LR_T2", "LR_T3"))
  meta_df$subject <- sub("_T[123]$", "", meta_df$Col_ID)

  design <- model.matrix(~ 0 + group, data = meta_df)
  colnames(design) <- gsub("^group", "", colnames(design))

  dupcor <- duplicateCorrelation(mat, design, block = meta_df$subject)
  fit <- lmFit(mat, design, block = meta_df$subject,
               correlation = dupcor$consensus.correlation)

  cm <- makeContrasts(
    Acute_LR = LR_T2 - LR_T1,
    levels = design)
  fit2 <- contrasts.fit(fit, cm)
  fit2 <- eBayes(fit2)

  topTable(fit2, coef = "Acute_LR", number = Inf, sort.by = "none")
}

# --- Run for each method ---
results <- list()

for (mname in names(imp_list)) {
  cat(sprintf("  Downstream: %s... ", mname))
  imp_mat <- imp_list[[mname]]

  # Ensure metadata alignment
  meta_df <- as.data.frame(meta)
  rownames(meta_df) <- meta_df$Col_ID

  # Run limma
  tt <- tryCatch(
    run_limma_benchmark(imp_mat, meta_df),
    error = function(e) { cat(sprintf("FAILED: %s\n", e$message)); NULL }
  )
  if (is.null(tt)) next

  # Gene names for matching
  tt$gene <- rownames(tt)

  # DEP count (FDR < 0.05)
  dep_count <- sum(tt$adj.P.Val < 0.05, na.rm = TRUE)

  # Store logFC for later rho computation
  results[[mname]] <- list(
    tt = tt,
    dep_count = dep_count
  )
  cat(sprintf("DEP=%d\n", dep_count))
}

# --- Compute FC rho and NES rho relative to Non_imputed ---
ref_name <- "Non_imputed"
if (!ref_name %in% names(results)) {
  stop("Non_imputed results not found — cannot compute relative metrics")
}
ref_tt <- results[[ref_name]]$tt
ref_logfc <- setNames(ref_tt$logFC, ref_tt$gene)

# Reference NES — filter out non-finite rank stats
ref_ranks <- setNames(-log10(ref_tt$P.Value) * sign(ref_tt$logFC), ref_tt$gene)
ref_ranks <- ref_ranks[is.finite(ref_ranks)]
ref_fgsea <- fgsea(pathways = gene_sets, stats = ref_ranks, minSize = 15, maxSize = 500)
ref_nes <- setNames(ref_fgsea$NES, ref_fgsea$pathway)

downstream_rows <- list()

for (mname in names(results)) {
  tt <- results[[mname]]$tt
  method_logfc <- setNames(tt$logFC, tt$gene)

  # FC Spearman rho (all shared genes)
  shared <- intersect(names(ref_logfc), names(method_logfc))
  fc_rho <- cor(ref_logfc[shared], method_logfc[shared], method = "spearman", use = "complete.obs")

  # NES rho — filter out non-finite rank stats
  method_ranks <- setNames(-log10(tt$P.Value) * sign(tt$logFC), tt$gene)
  method_ranks <- method_ranks[is.finite(method_ranks)]
  method_fgsea <- fgsea(pathways = gene_sets, stats = method_ranks, minSize = 15, maxSize = 500)
  method_nes <- setNames(method_fgsea$NES, method_fgsea$pathway)
  shared_pw <- intersect(names(ref_nes), names(method_nes))
  nes_rho <- cor(ref_nes[shared_pw], method_nes[shared_pw], method = "spearman", use = "complete.obs")

  downstream_rows[[mname]] <- data.frame(
    method = mname,
    fc_rho = fc_rho,
    nes_rho = nes_rho,
    dep_count = results[[mname]]$dep_count,
    stringsAsFactors = FALSE
  )
  cat(sprintf("  %s: FC_rho=%.3f NES_rho=%.3f DEP=%d\n",
              mname, fc_rho, nes_rho, results[[mname]]$dep_count))
}

downstream_df <- do.call(rbind, downstream_rows)
rownames(downstream_df) <- NULL
write.csv(downstream_df, file.path(BENCH_DIR, "02_downstream.csv"), row.names = FALSE)
cat(sprintf("Wrote %d rows to 02_downstream.csv\n", nrow(downstream_df)))

# HRvLR WGCNA Runner — produces network and module artifacts for F06/F08/F09
# Inputs:  02_Normalization/imputation/c_data/DAList_imputed_missforest.rds
# Outputs: 04_Figures/F06/c_data/wgcna/  (network, sft_summary, stratified trait cors)
#          04_Figures/F06/c_data/        (panel-ready: MEs, kME, module_colors, datExpr,
#                                         module_assignments, mod_bio_labels)
#
# Mirrors YvO_WGCNA_run.R parameters (Cahill 2018 / Langfelder & Horvath 2008):
#   networkType = "signed", TOMType = "signed",
#   minModuleSize = 30, mergeCutHeight = 0.25,
#   pickSoftThreshold over powers 1:20, R^2 > 0.87 cutoff (default 6 fallback).
# Stratification: HR and LR responder groups (mirrors YvO's Young/Old split).

suppressPackageStartupMessages({
  library(WGCNA)
  library(tidyverse)
})

allowWGCNAThreads()
set.seed(42)
setwd(rprojroot::find_rstudio_root_file())

DALIST_RDS <- "02_Normalization/imputation/c_data/DAList_imputed_missforest.rds"
PANEL_DIR <- "04_Figures/F06/c_data"
DATA_DIR <- file.path(PANEL_DIR, "wgcna")
REPORT_DIR <- "04_Figures/F06/b_reports/supp/01_QC"

dir.create(PANEL_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(REPORT_DIR, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(DALIST_RDS))

dal <- readRDS(DALIST_RDS)
df <- dplyr::bind_cols(
  tibble::as_tibble(dal$annotation[, c("uniprot_id", "protein", "gene", "description")]),
  tibble::as_tibble(dal$data)
)
ann_cols <- c("uniprot_id", "protein", "gene", "description")
ann <- df[, ann_cols]
samp_names <- setdiff(names(df), ann_cols)
mat <- as.matrix(df[, samp_names])
rownames(mat) <- ann$uniprot_id

datExpr <- t(mat)

dal_meta <- as.data.frame(dal$metadata)

meta <- tibble(
  sample_id  = dal_meta$Col_ID,
  subject    = dal_meta$Subject_ID, # HRvLR IDs are group-prefixed (HR_S03 / LR_S04)
  group      = dal_meta$Group,
  timepoint  = dal_meta$Timepoint,
  time       = dal_meta$Timepoint, # alias for downstream supp panels
  group_time = dal_meta$Group_Time
)
# Pass phenotype columns through from dal_meta for downstream supp panels
pheno_cols <- intersect(
  c(
    "X1RM_Leg_Pre", "mCSA_Pre", "fCSA_Mixed_Pre", "fCSA_Type_I_Pre",
    "fCSA_Type_II_Pre", "MyoVision_fCSA_mixed_Pre", "MyoVision_fCSA_Type_I__Pre",
    "ACCUM_VL", "X1RM._Ext_Pre", "COMP.HYPERTROPHY"
  ),
  names(dal_meta)
)
for (pc in pheno_cols) {
  meta[[pc]] <- dal_meta[[pc]][match(meta$sample_id, dal_meta$Col_ID)]
}

gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  ann <- ann |> filter(uniprot_id %in% colnames(datExpr))
  message(sprintf(
    "After goodSamplesGenes: %d samples x %d proteins",
    nrow(datExpr), ncol(datExpr)
  ))
}

cor <- WGCNA::cor

powers <- 1:20
sft <- pickSoftThreshold(datExpr,
  powerVector = powers,
  networkType = "signed", verbose = 2
)
saveRDS(sft$fitIndices, file.path(DATA_DIR, "sft_fitIndices.rds"))

r2_values <- -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq
power_idx <- which(r2_values > 0.87)[1]
soft_power <- if (!is.na(power_idx)) powers[power_idx] else 6L
message(sprintf("Soft power: %d (R^2 = %.3f)", soft_power, r2_values[soft_power]))

png(file.path(REPORT_DIR, "SUPP_soft_threshold.png"),
  width = 3000, height = 1500, res = 300
)
par(mfrow = c(1, 2))
plot(sft$fitIndices$Power, r2_values,
  xlab = "Soft Threshold (power)",
  ylab = expression(paste("Scale Free Topology Model Fit (", R^2, ")")),
  main = "Scale independence", type = "n"
)
text(sft$fitIndices$Power, r2_values, labels = powers, cex = 0.9, col = "red")
abline(h = 0.85, col = "red", lty = 2)
plot(sft$fitIndices$Power, sft$fitIndices$mean.k.,
  xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
  main = "Mean connectivity", type = "n"
)
text(sft$fitIndices$Power, sft$fitIndices$mean.k.,
  labels = powers,
  cex = 0.9, col = "red"
)
dev.off()

net <- blockwiseModules(
  datExpr,
  power             = soft_power,
  networkType       = "signed",
  TOMType           = "signed",
  minModuleSize     = 30,
  mergeCutHeight    = 0.25,
  numericLabels     = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs          = FALSE,
  verbose           = 3
)

module_colors <- labels2colors(net$colors)
n_modules <- length(unique(net$colors)) - (0 %in% net$colors)
message(sprintf("Modules detected: %d (+ grey/unassigned)", n_modules))

MEs <- moduleEigengenes(datExpr, colors = module_colors)$eigengenes
MEs <- orderMEs(MEs)
kME <- signedKME(datExpr, MEs)

module_df <- tibble(
  uniprot_id   = colnames(datExpr),
  module_color = module_colors,
  module_num   = net$colors
) |>
  left_join(ann |> dplyr::select(uniprot_id, gene), by = "uniprot_id")

mod_sizes <- sort(table(module_colors[module_colors != "grey"]),
  decreasing = TRUE
)
mod_bio_labels <- tibble(
  module_color = names(mod_sizes),
  module_id = paste0("M", seq_along(mod_sizes)),
  bio_label = NA_character_, # hand-curated downstream
  n_proteins = as.integer(mod_sizes),
  display_label = paste0(
    "M", seq_along(mod_sizes), ": ",
    names(mod_sizes), " (n=", as.integer(mod_sizes), ")"
  )
)

# Stratified trait correlations (HR vs LR), mirrors YvO Young/Old split.
# Traits drawn from dal_meta numeric columns (HRvLR phenotype columns are stored
# inline in the DAList metadata — _Pre suffix marks T1 baseline measurements).
trait_cols <- names(dal_meta)[vapply(dal_meta, is.numeric, logical(1))]
trait_cols <- setdiff(trait_cols, c("Subject_ID", "subject", "sample_id"))

cor <- stats::cor

# Baseline (T1) per stratum
write_baseline_cors <- function(stratum_label, stratum_subjects) {
  sub_idx <- meta$subject %in% stratum_subjects & meta$timepoint == "T1"
  if (sum(sub_idx) < 5 || length(trait_cols) == 0) {
    return(invisible(NULL))
  }
  me_sub <- MEs[meta$sample_id[sub_idx], , drop = FALSE]
  trait_mat <- as.data.frame(dal_meta[match(
    meta$sample_id[sub_idx],
    dal_meta$Col_ID
  ), trait_cols, drop = FALSE])
  cor_mat <- cor(me_sub, trait_mat, use = "pairwise.complete.obs")
  pval_mat <- corPvalueStudent(cor_mat, nrow(me_sub))
  bh_mat <- apply(pval_mat, 2, p.adjust, method = "BH")
  to_df <- function(m) as.data.frame(m) |> rownames_to_column("module")
  write_csv(
    to_df(cor_mat),
    file.path(
      DATA_DIR,
      sprintf("wgcna_baseline_trait_correlations_%s.csv", stratum_label)
    )
  )
  write_csv(
    to_df(pval_mat),
    file.path(
      DATA_DIR,
      sprintf("wgcna_baseline_trait_pvalues_raw_%s.csv", stratum_label)
    )
  )
  write_csv(
    to_df(bh_mat),
    file.path(
      DATA_DIR,
      sprintf("wgcna_baseline_trait_pvalues_bh_%s.csv", stratum_label)
    )
  )
}

# Change (T2 - T1 paired) per stratum: ME deltas vs trait deltas where available
write_change_cors <- function(stratum_label, stratum_subjects) {
  sub_t1 <- meta$subject %in% stratum_subjects & meta$timepoint == "T1"
  sub_t2 <- meta$subject %in% stratum_subjects & meta$timepoint == "T2"
  if (sum(sub_t1) < 3 || sum(sub_t2) < 3 || length(trait_cols) == 0) {
    return(invisible(NULL))
  }
  me_t1 <- MEs[meta$sample_id[sub_t1], , drop = FALSE]
  me_t2 <- MEs[meta$sample_id[sub_t2], , drop = FALSE]
  rownames(me_t1) <- meta$subject[sub_t1]
  rownames(me_t2) <- meta$subject[sub_t2]
  common <- intersect(rownames(me_t1), rownames(me_t2))
  if (length(common) < 3) {
    return(invisible(NULL))
  }
  delta_me <- me_t2[common, , drop = FALSE] - me_t1[common, , drop = FALSE]
  trait_t1 <- as.data.frame(dal_meta[match(paste0(common, "_T1"), dal_meta$Col_ID),
    trait_cols,
    drop = FALSE
  ])
  trait_t2 <- as.data.frame(dal_meta[match(paste0(common, "_T2"), dal_meta$Col_ID),
    trait_cols,
    drop = FALSE
  ])
  delta_trait <- as.matrix(trait_t2) - as.matrix(trait_t1)
  cor_mat <- cor(delta_me, delta_trait, use = "pairwise.complete.obs")
  pval_mat <- corPvalueStudent(cor_mat, length(common))
  bh_mat <- apply(pval_mat, 2, p.adjust, method = "BH")
  to_df <- function(m) as.data.frame(m) |> rownames_to_column("module")
  write_csv(
    to_df(cor_mat),
    file.path(
      DATA_DIR,
      sprintf("wgcna_change_trait_correlations_%s.csv", stratum_label)
    )
  )
  write_csv(
    to_df(pval_mat),
    file.path(
      DATA_DIR,
      sprintf("wgcna_change_trait_pvalues_raw_%s.csv", stratum_label)
    )
  )
  write_csv(
    to_df(bh_mat),
    file.path(
      DATA_DIR,
      sprintf("wgcna_change_trait_pvalues_bh_%s.csv", stratum_label)
    )
  )
}

write_baseline_cors("hr", meta$subject[meta$group == "HR"])
write_baseline_cors("lr", meta$subject[meta$group == "LR"])
write_change_cors("hr", meta$subject[meta$group == "HR"])
write_change_cors("lr", meta$subject[meta$group == "LR"])

# --- LMM contrasts: ME ~ group_time + (1|subject), KR df ---
# Contrasts mirror YvO's pattern: per-stratum training + acute effects.
# 6-level group_time: HR_T1, HR_T2, HR_T3, LR_T1, LR_T2, LR_T3
lmm_contrasts <- list(
  Training_HR = c(-1, 1, 0, 0, 0, 0),
  Training_LR = c(0, 0, 0, -1, 1, 0),
  Acute_HR    = c(-1, 0, 1, 0, 0, 0),
  Acute_LR    = c(0, 0, 0, -1, 0, 1)
)
gt_levels <- c("HR_T1", "HR_T2", "HR_T3", "LR_T1", "LR_T2", "LR_T3")
gt_levels <- intersect(gt_levels, unique(meta$group_time))

if (length(gt_levels) >= 4) {
  lmm_data <- meta |>
    mutate(group_time = factor(group_time, levels = gt_levels))
  lmm_rows <- list()
  for (mod in colnames(MEs)) {
    lmm_data[[mod]] <- MEs[lmm_data$sample_id, mod]
    fit <- tryCatch(suppressWarnings(
      lme4::lmer(as.formula(paste0("`", mod, "` ~ group_time + (1 | subject)")),
        data = lmm_data, REML = TRUE
      )
    ), error = function(e) NULL)
    if (is.null(fit)) next
    emm <- emmeans::emmeans(fit, ~group_time)
    for (cname in names(lmm_contrasts)) {
      cv <- lmm_contrasts[[cname]][seq_along(gt_levels)]
      ctr <- emmeans::contrast(emm, list(ctr = cv))
      s <- summary(ctr, ddf = "Kenward-Roger")
      r_eq <- sign(s$estimate) * sqrt(s$t.ratio^2 / (s$t.ratio^2 + s$df))
      lmm_rows <- c(lmm_rows, list(tibble(
        module   = mod,
        contrast = cname,
        estimate = round(s$estimate, 5),
        SE       = round(s$SE, 5),
        df       = round(s$df, 2),
        t_ratio  = round(s$t.ratio, 4),
        p_raw    = s$p.value,
        r_equiv  = round(r_eq, 4)
      )))
    }
  }
  lmm_df <- bind_rows(lmm_rows)
  lmm_df$p_bh <- p.adjust(lmm_df$p_raw, method = "BH")
  write_csv(lmm_df, file.path(DATA_DIR, "wgcna_lmm_contrast_audit.csv"))
  message(sprintf(
    "LMM audit: %d tests (%d modules x %d contrasts)",
    nrow(lmm_df), length(colnames(MEs)), length(lmm_contrasts)
  ))
}

sft_summary <- tibble(
  selected_power    = soft_power,
  R_squared         = r2_values[soft_power],
  mean_connectivity = sft$fitIndices$mean.k.[soft_power],
  n_proteins        = ncol(datExpr),
  n_samples         = nrow(datExpr)
)

# Per-module ORA enrichment (consumed by F09/panel_B_triptych.R)
source("04_Figures/shared/pathway_utils.R")
source("04_Figures/shared/style.R")
bg_genes <- unique(ann$gene[!is.na(ann$gene) & nzchar(ann$gene)])
pw_collection <- build_pathway_collection(min_size = 15, max_size = 500)
enrich_rows <- list()
for (mc in setdiff(unique(module_colors), "grey")) {
  mod_genes <- unique(module_df$gene[module_df$module_color == mc])
  mod_genes <- mod_genes[!is.na(mod_genes) & nzchar(mod_genes)]
  if (length(mod_genes) < 5) next
  res <- tryCatch(
    run_ora_deduplicated(
      genes = mod_genes, universe = bg_genes,
      pathways = pw_collection, jaccard_cutoff = 0.5,
      min_size = 15, max_size = 500, padj_cutoff = 0.10
    ),
    error = function(e) NULL
  )
  if (is.null(res) || nrow(res) == 0) next
  res$module <- mc
  res$Description <- clean_pathway_name(res$pathway)
  res$ID <- res$pathway
  res$p.adjust <- res$padj
  res$Count <- res$overlap
  res$geneID <- vapply(
    res$overlapGenes,
    function(g) paste(g, collapse = "/"), character(1)
  )
  res$overlapGenes <- NULL
  enrich_rows <- c(enrich_rows, list(res))
}
enrich_df <- bind_rows(enrich_rows)
write_csv(enrich_df, file.path(DATA_DIR, "wgcna_module_enrichment.csv"))
message(sprintf(
  "Module enrichment: %d pathways across %d modules",
  nrow(enrich_df), length(setdiff(unique(module_colors), "grey"))
))

# WGCNA network artifacts (DATA_DIR)
saveRDS(net, file.path(DATA_DIR, "wgcna_network.rds"))
write_csv(sft_summary, file.path(DATA_DIR, "wgcna_sft_summary.csv"))
write_csv(module_df, file.path(DATA_DIR, "wgcna_module_assignments.csv"))

# Panel-ready artifacts (PANEL_DIR)
imp_mat <- t(datExpr) # proteins x samples for downstream triptych panels
saveRDS(MEs, file.path(PANEL_DIR, "MEs.rds"))
saveRDS(kME, file.path(PANEL_DIR, "kME_all.rds"))
saveRDS(datExpr, file.path(PANEL_DIR, "datExpr.rds"))
saveRDS(imp_mat, file.path(PANEL_DIR, "imp_mat.rds"))
saveRDS(module_colors, file.path(PANEL_DIR, "module_colors.rds"))
write_csv(module_df, file.path(PANEL_DIR, "wgcna_module_assignments.csv"))
write_csv(mod_bio_labels, file.path(PANEL_DIR, "mod_bio_labels.csv"))
write_csv(meta, file.path(PANEL_DIR, "meta.csv"))
write_csv(ann, file.path(PANEL_DIR, "imp_annotations.csv"))

key_modules <- head(mod_bio_labels$module_color, 5)
writeLines(key_modules, file.path(PANEL_DIR, "key_modules.txt"))

message(sprintf(
  "Done: %d modules x %d proteins x %d samples (power = %d)",
  n_modules, ncol(datExpr), nrow(datExpr), soft_power
))

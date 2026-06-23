# Paired-design WGCNA on the pi-gated proteome (298 proteins × 45 samples).
# Within-subject centering removes between-subject identity from modules
# before network construction. See Li et al. 2018 Sci Rep for rationale.

pacman::p_load(here, WGCNA)
# Load data — WGCNA masks stats::cor, so qualify downstream if needed
source(here("05_Clustering", "_shared_inputs.R"))
i <- load_clustering_inputs()

# within-subject centering
center_within_subject <- function(abund, meta) {
  subjects <- unique(meta$subject)
  out <- abund
  for (s in subjects) {
    cols <- meta$sample_id[meta$subject == s]
    sub_mat <- abund[, cols, drop = FALSE]
    out[, cols] <- sub_mat - rowMeans(sub_mat)
  }
  out
}

abund_c <- center_within_subject(i$abund, i$meta)
expr <- t(abund_c) # 45 samples × 298 proteins

# soft threshold
powers <- 1:20
sft <- pickSoftThreshold(
  expr,
  powerVector = powers,
  networkType = "signed",
  verbose = 0
)

fit <- sft$fitIndices
r2 <- -sign(fit$slope) * fit$SFT.R.sq

candidates <- which(r2 >= 0.8)
if (length(candidates)) {
  chosen_power <- powers[candidates[1]]
  fallback <- FALSE
} else {
  # fallback: scale-free R² never reached 0.8; use signed-network default
  chosen_power <- 12L
  fallback <- TRUE
  message(
    "No power reached scale-free R² ≥ 0.8; ",
    "using signed-network fallback power = ", chosen_power
  )
}

message(
  "Chosen soft-threshold power: ", chosen_power,
  " (fallback: ", fallback, ")"
)

pdf(here("05_Clustering", "a_wgcna_paired", "b_reports", "soft_threshold.pdf"),
  width = 9, height = 4.5
)
par(mfrow = c(1, 2))
plot(fit$Power, r2,
  xlab = "Soft Threshold (power)", ylab = "Scale Free Topology R²",
  type = "b", main = "Scale-free fit", pch = 16
)
abline(h = 0.8, col = "red", lty = 2)
abline(v = chosen_power, col = "steelblue", lty = 2)
plot(fit$Power, fit$mean.k.,
  xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
  type = "b", main = "Mean connectivity", pch = 16
)
abline(v = chosen_power, col = "steelblue", lty = 2)
dev.off()

# blockwise modules
net <- blockwiseModules(
  expr,
  networkType = "signed",
  power = chosen_power,
  minModuleSize = 20,
  mergeCutHeight = 0.25,
  deepSplit = 2,
  verbose = 0
)

pdf(here("05_Clustering", "a_wgcna_paired", "b_reports", "dendro_modules.pdf"),
  width = 10, height = 5
)
plotDendroAndColors(
  net$dendrograms[[1]],
  net$colors,
  "Module",
  dendroLabels = FALSE,
  hang = 0.03,
  main = "Gene dendrogram and module colors"
)
dev.off()

# eigengenes + sign fix
mes_all <- moduleEigengenes(expr, colors = net$colors)$eigengenes
stopifnot(rownames(mes_all) == i$meta$sample_id)

# drop the grey (unassigned) eigengene from the ME table
grey_col <- grep("^MEgrey$", names(mes_all), value = TRUE)
mes_real <- mes_all[, setdiff(names(mes_all), grey_col), drop = FALSE]

# sign-fix: positive ME → higher mean centered abundance for that module
for (mod_col in names(mes_real)) {
  mod_label <- sub("^ME", "", mod_col)
  in_module <- net$colors == mod_label
  mod_mean <- colMeans(abund_c[in_module, , drop = FALSE])
  if (stats::cor(mes_real[[mod_col]], mod_mean) < 0) {
    mes_real[[mod_col]] <- -mes_real[[mod_col]]
  }
}

module_labels <- sub("^ME", "", names(mes_real))
n_modules <- length(module_labels)
message("Non-grey modules: ", n_modules)
message("Grey (unassigned) proteins: ", sum(net$colors == "grey"))

message("\n-- Sensitivity: powers near chosen --")
for (p in c(chosen_power - 2L, chosen_power, chosen_power + 2L)) {
  if (p < 1 || p > 20) next
  net_s <- blockwiseModules(
    expr,
    networkType    = "signed",
    power          = p,
    minModuleSize  = 20,
    mergeCutHeight = 0.25,
    deepSplit      = 2,
    verbose        = 0
  )
  message(
    "  power=", p, " → modules=",
    length(setdiff(unique(net_s$colors), "grey"))
  )
}

message("\n-- Sensitivity: mergeCutHeight --")
for (mch in c(0.15, 0.25, 0.35)) {
  net_s <- blockwiseModules(
    expr,
    networkType    = "signed",
    power          = chosen_power,
    minModuleSize  = 20,
    mergeCutHeight = mch,
    deepSplit      = 2,
    verbose        = 0
  )
  message(
    "  mergeCutHeight=", mch, " → modules=",
    length(setdiff(unique(net_s$colors), "grey"))
  )
}

# write artifacts
membership <- data.frame(
  protein_id = colnames(expr),
  group_id = net$colors,
  membership_weight = 1L,
  stringsAsFactors = FALSE
)

eigengene_rows <- lapply(seq_len(n_modules), function(j) {
  data.frame(
    sample_id = i$meta$sample_id,
    subject = i$meta$subject,
    group_arm = i$meta$group_arm,
    timepoint = i$meta$timepoint,
    group_id = module_labels[j],
    ME = mes_real[[j]],
    stringsAsFactors = FALSE
  )
})
eigengene <- do.call(rbind, eigengene_rows)

write.csv(
  membership,
  here("05_Clustering", "a_wgcna_paired", "c_data", "membership.csv"),
  row.names = FALSE
)
write.csv(
  eigengene,
  here("05_Clustering", "a_wgcna_paired", "c_data", "eigengene.csv"),
  row.names = FALSE
)

message(
  "Done. membership.csv: ", nrow(membership), " rows | eigengene.csv: ",
  nrow(eigengene), " rows (", n_modules, " modules × 45 samples)"
)

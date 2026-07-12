pacman::p_load(WGCNA)

# abund: proteins x samples. Returns the single-source module fit.
fit_modules <- function(abund, sample_meta) {
  WGCNA::disableWGCNAThreads()
  expr <- t(abund)
  gsg <- WGCNA::goodSamplesGenes(expr, verbose = 0)
  if (!gsg$allOK) expr <- expr[gsg$goodSamples, gsg$goodGenes]

  powers <- 1:20
  sft <- WGCNA::pickSoftThreshold(expr, powerVector = powers, networkType = "signed", verbose = 0)
  r2 <- -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq
  hit <- which(r2 > 0.87)
  chosen_power <- if (length(hit)) powers[hit[1]] else 12L

  set.seed(42)
  net <- WGCNA::blockwiseModules(
    expr,
    networkType = "signed", TOMType = "signed", power = chosen_power,
    minModuleSize = 30, mergeCutHeight = 0.25, numericLabels = TRUE,
    pamRespectsDendro = FALSE, saveTOMs = FALSE, randomSeed = 42, verbose = 0
  )
  colors <- WGCNA::labels2colors(net$colors)
  names(colors) <- colnames(expr)

  mes <- WGCNA::moduleEigengenes(expr, colors = colors)$eigengenes
  mes <- mes[, setdiff(colnames(mes), "MEgrey"), drop = FALSE]
  rownames(mes) <- rownames(expr)

  membership <- data.frame(protein_id = colnames(expr), module = colors, stringsAsFactors = FALSE)
  list(colors = colors, eigengenes = mes, membership = membership, chosen_power = chosen_power)
}

# Co-expression modules: the module-by-sample eigengene matrix and each protein's membership,
# which every later network step reads.
#
# WGCNA treats samples as independent, but these 45 samples are 16 subjects measured up to three
# times, and subject identity drives the leading components of raw abundance. Modules are
# defined on subject-centred abundance and scored on raw abundance, so between-arm differences
# stay testable. A subject with one biopsy centres to a column of zeros and carries no
# within-subject information, so it is left out of the definition and kept in the scoring.
# Construction never sees a label or a phenotype. WGCNA needs a complete matrix, so this reads the
# imputed one.
#
# Parameters follow the WGCNA FAQ's two recommended departures from default, a signed network and
# biweight midcorrelation. minModuleSize is 30 against a default of min(20, n/2), and maxBlockSize
# 2500 keeps all 1,900 proteins in one block; the rest are package defaults.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(WGCNA)
})

stage <- here("04_Network", "01_Build_Modules")
out <- file.path(stage, "c_data")
figure_dir <- file.path(stage, "b_reports")
for (path in c(out, figure_dir)) dir.create(path, recursive = TRUE, showWarnings = FALSE)

inputs <- c(imputed = "01_Preprocess/03_Imputation/c_data/DAList_imputed.rds")
paths <- map_chr(inputs, here)
if (!all(file.exists(paths))) {
  stop("Run 01_Preprocess first. Missing: ", paste(inputs[!file.exists(paths)], collapse = ", "))
}
imputed <- readRDS(paths[["imputed"]])

abund <- as.matrix(imputed$data)
meta <- as_tibble(imputed$metadata) |>
  select(sample_id, subject, arm, timepoint)
stopifnot(identical(colnames(abund), meta$sample_id))
disableWGCNAThreads()

# maxPOutliers = 0.05 is the FAQ's "strongly recommend"; the default of 1 disables the outlier
# guard. pearsonFallback covers proteins with zero MAD, which imputation can produce.
bicor_args <- list(corType = "bicor", maxPOutliers = 0.05, pearsonFallback = "individual")
repeated <- meta |>
  add_count(subject) |>
  filter(n > 1)
centred <- t(apply(abund[, repeated$sample_id], 1, \(x) x - ave(x, repeated$subject)))
expr <- t(centred)

# The power is pickSoftThreshold's estimate: the lowest whose scale-free fit R2 clears 0.85.
# The WGCNA FAQ's signed table (under 20 samples 18, 20-30 16, 31-40 14, over 40 12) is the
# fallback when no power clears it.
sft <- pickSoftThreshold(expr,
  powerVector = 1:20, networkType = "signed",
  corFnc = "bicor", corOptions = list(maxPOutliers = 0.05), verbose = 0
)
faq_power <- c(18L, 16L, 14L, 12L)[findInterval(nrow(expr), c(20, 31, 41)) + 1]
power <- coalesce(sft$powerEstimate, faq_power)
fit_indices <- as_tibble(sft$fitIndices) |>
  mutate(signed_r2 = -sign(slope) * SFT.R.sq, chosen = Power == power)

# randomSeed seeds the clustering; a set.seed() here would be overridden.
net <- do.call(blockwiseModules, c(
  list(expr,
    power = power, networkType = "signed", TOMType = "signed",
    deepSplit = 2, minModuleSize = 30, mergeCutHeight = 0.15,
    pamStage = TRUE, pamRespectsDendro = TRUE, maxBlockSize = 2500,
    numericLabels = FALSE, randomSeed = 54321, verbose = 0
  ),
  bicor_args
))
colours <- set_names(net$colors, colnames(expr))

# moduleEigengenes aligns each eigengene with its module's average expression by default.
me <- moduleEigengenes(t(abund), colors = colours, excludeGrey = TRUE)$eigengenes
rownames(me) <- colnames(abund)
eigengenes <- t(as.matrix(me))
rownames(eigengenes) <- sub("^ME", "", rownames(eigengenes))

kme <- signedKME(t(abund), me,
  corFnc = "bicor", corOptions = "maxPOutliers = 0.05, pearsonFallback = 'individual'"
)
colnames(kme) <- sub("^kME", "", colnames(kme))
gene_of <- set_names(imputed$annotation$gene, imputed$annotation$uniprot_id)
membership <- tibble(
  uniprot_id = names(colours), gene = unname(gene_of[names(colours)]), module = unname(colours)
) |>
  mutate(kme = map2_dbl(uniprot_id, module, \(id, m) if (m == "grey") NA_real_ else kme[id, m]))

eigengene_long <- as_tibble(t(eigengenes), rownames = "sample_id") |>
  pivot_longer(-sample_id, names_to = "module", values_to = "eigengene") |>
  left_join(meta, by = "sample_id")

# Whether subject identity still drives a scored eigengene. Raw R2 is inflated by 15 subject
# degrees of freedom on 45 samples, so the intraclass correlation is the number reported.
subject_icc <- eigengene_long |>
  summarise(
    icc = {
      v <- as.data.frame(lme4::VarCorr(suppressMessages(
        lme4::lmer(eigengene ~ 1 + (1 | subject))
      )))$vcov
      v[1] / (v[1] + v[2])
    },
    .by = module
  )
module_summary <- count(filter(membership, module != "grey"), module, name = "n_proteins") |>
  left_join(subject_icc, by = "module") |>
  arrange(desc(n_proteins))
print(as.data.frame(module_summary), digits = 3)
message(sprintf(
  "power %d (signed R2 %.3f, mean k %.1f); %d modules, %d of %d proteins unassigned; %s",
  power, fit_indices$signed_r2[power], fit_indices$mean.k.[power], nrow(eigengenes),
  sum(colours == "grey"), length(colours),
  sprintf("%d of %d samples define the modules", nrow(expr), ncol(abund))
))


# ---- figures -------------------------------------------------------------------------------

caption_theme <- theme(plot.caption = element_text(size = 7, colour = "grey45", hjust = 0))

threshold_figure <- fit_indices |>
  select(Power, `signed scale-free R2` = signed_r2, `mean connectivity` = mean.k.) |>
  pivot_longer(-Power) |>
  ggplot(aes(Power, value)) +
  geom_line(colour = "grey60") +
  geom_point(aes(colour = Power == power), size = 2) +
  geom_hline(
    data = tibble(name = "signed scale-free R2", y = 0.85), aes(yintercept = y),
    linetype = "dashed"
  ) +
  facet_wrap(~name, scales = "free_y") +
  scale_colour_manual(values = c(`FALSE` = "grey30", `TRUE` = "#B2182B"), guide = "none") +
  labs(
    x = "soft power", y = NULL, title = "Soft-threshold choice",
    subtitle = sprintf("WGCNA signed network, bicor, subject-centred matrix; power %d", power),
    caption = paste(
      "Left: mean connectivity at each power. Right: signed scale-free fit; the dashed line is",
      "the 0.85 criterion and the red point the chosen power.",
      "Table: c_data/01_build_modules.xlsx, soft_threshold."
    )
  ) +
  theme_minimal(base_size = 10) +
  caption_theme

module_colours <- set_names(module_summary$module, module_summary$module)
size_figure <- module_summary |>
  mutate(module = factor(module, levels = rev(module))) |>
  ggplot(aes(n_proteins, module, fill = module)) +
  geom_col(colour = "grey30", linewidth = 0.2) +
  geom_text(aes(label = sprintf("ICC %.2f", icc)), hjust = -0.15, size = 3) +
  scale_fill_manual(values = module_colours, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(
    x = "proteins", y = NULL, title = "Module sizes and subject dependence",
    subtitle = sprintf(
      "%d modules; %d of %d proteins unassigned", nrow(module_summary),
      sum(colours == "grey"), length(colours)
    ),
    caption = paste(
      "Bars: proteins per module. Label: intraclass correlation of the eigengene across a",
      "subject's biopsies (lme4, 1 | subject): the share of eigengene variance between subjects.",
      "Table: c_data/01_build_modules.xlsx, module_summary."
    )
  ) +
  theme_minimal(base_size = 10) +
  caption_theme

trajectory <- eigengene_long |>
  summarise(
    mean = mean(eigengene), se = sd(eigengene) / sqrt(n()),
    .by = c(module, arm, timepoint)
  )
trajectory_figure <- ggplot(trajectory, aes(timepoint, mean, colour = arm, group = arm)) +
  geom_line(
    data = eigengene_long, aes(y = eigengene, group = subject),
    alpha = 0.2, linewidth = 0.3
  ) +
  geom_line(linewidth = 0.8) +
  geom_pointrange(aes(ymin = mean - se, ymax = mean + se), size = 0.25) +
  facet_wrap(~ factor(module, levels = module_summary$module), nrow = 3) +
  scale_colour_manual(values = c(HR = "#2166AC", LR = "#B2182B"), name = NULL) +
  labs(
    x = NULL, y = "eigengene", title = "Eigengene trajectories",
    subtitle = "module eigengene per biopsy, scored on raw abundance",
    caption = paste(
      "Thin lines: one subject's biopsies. Thick line and range: arm mean and standard error",
      "at each timepoint. Table: c_data/01_build_modules.xlsx, eigengenes."
    )
  ) +
  theme_minimal(base_size = 10) +
  caption_theme

pdf(file.path(figure_dir, "01_build_modules_figures.pdf"), width = 11, height = 8.5)
walk(list(threshold_figure, size_figure, trajectory_figure), print)
invisible(dev.off())

saveRDS(
  list(
    eigengenes = eigengenes, membership = membership, meta = meta, module_summary = module_summary
  ),
  file.path(out, "modules.rds"),
  compress = "xz"
)
sheets <- list(
  module_summary = module_summary, membership = membership,
  eigengenes = as_tibble(t(eigengenes), rownames = "sample_id"),
  soft_threshold = fit_indices,
  input_manifest = tibble(
    input = names(inputs), path = unname(inputs), md5 = unname(tools::md5sum(paths))
  ),
  package_versions = sessioninfo::package_info("loaded", dependencies = FALSE) |>
    as_tibble() |>
    select(package, version = loadedversion, source)
)
read_me <- tibble(sheet = names(sheets), holds = c(
  "Proteins per module and the subject ICC of its eigengene.",
  "Every protein's module and its kME in that module.",
  "Module eigengene per sample, scored on raw abundance.",
  "pickSoftThreshold fit indices per power; chosen marks the power used.",
  "Files read, with md5.",
  "Packages loaded at run time."
))
writexl::write_xlsx(c(list(read_me = read_me), sheets), file.path(out, "01_build_modules.xlsx"))
message("wrote modules.rds, 01_build_modules.xlsx and 3 figures")

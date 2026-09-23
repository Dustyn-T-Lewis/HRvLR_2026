# Co-expression modules: the module-by-sample eigengene matrix and each protein's membership,
# which every later network step reads.
#
# WGCNA correlates proteins across samples as if they were independent. These 45 samples are 16
# subjects measured up to three times, and on raw abundance subject identity drives the leading
# components, so modules built there would encode who a biopsy came from. Modules are therefore
# defined on abundance centred within subject, which leaves how proteins move together inside a
# person, and eigengenes are scored on raw abundance, so between-arm differences survive to be
# tested. Construction never sees a label or a phenotype. WGCNA needs a complete matrix, so this
# reads the imputed one.
#
# Parameters follow the WGCNA FAQ's two recommended departures from default, a signed network and
# biweight midcorrelation, and are otherwise package defaults.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(WGCNA)
})

stage <- here("04_Network", "01_build_modules")
out <- file.path(stage, "c_data")
figure_dir <- file.path(stage, "b_reports")
for (path in c(out, figure_dir)) dir.create(path, recursive = TRUE, showWarnings = FALSE)

inputs <- c(imputed = "01_Preprocess/03_Imputation/c_data/DAList_imputed.rds")
paths <- map_chr(inputs, here)
if (!all(file.exists(paths))) {
  stop("Run 01_Preprocess first. Missing: ", paste(inputs[!file.exists(paths)], collapse = ", "))
}
imputed <- readRDS(paths[["imputed"]])
manifest <- tibble(
  input = names(inputs), path = unname(inputs), md5 = unname(tools::md5sum(paths))
)

abund <- as.matrix(imputed$data)
meta <- as_tibble(imputed$metadata) |>
  select(sample_id = Col_ID, subject = Subject_ID, arm = Group, timepoint = Timepoint)
stopifnot(identical(colnames(abund), meta$sample_id))
disableWGCNAThreads()

# maxPOutliers = 0.05 is the FAQ's "strongly recommend"; the default of 1 disables the outlier
# guard. pearsonFallback covers proteins with zero MAD, which imputation can produce.
bicor_args <- list(corType = "bicor", maxPOutliers = 0.05, pearsonFallback = "individual")
centred <- abund - t(apply(t(abund), 2, \(x) ave(x, meta$subject)))
expr <- t(centred)

# The power is the lowest whose signed scale-free fit clears pickSoftThreshold's default 0.85.
# The FAQ's sample-size table is the fallback when no power clears it.
sft <- pickSoftThreshold(expr,
  powerVector = 1:20, networkType = "signed",
  corFnc = "bicor", corOptions = list(maxPOutliers = 0.05), verbose = 0
)
power <- coalesce(sft$powerEstimate, if (nrow(expr) > 40) 12L else 16L)
fit_indices <- as_tibble(sft$fitIndices) |>
  mutate(signed_r2 = -sign(slope) * SFT.R.sq)

set.seed(42)
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
membership <- tibble(
  uniprot_id = names(colours), gene = imputed$annotation$gene[match(names(colours),
    imputed$annotation$uniprot_id)],
  module = unname(colours)
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
module_summary <- count(filter(membership, module != "grey"), module, name = "proteins") |>
  left_join(subject_icc, by = "module") |>
  arrange(desc(proteins))
print(as.data.frame(module_summary), digits = 3)
message(sprintf(
  "power %d (signed R2 %.3f, mean k %.1f); %d modules, %d of %d proteins unassigned",
  power, fit_indices$signed_r2[power], fit_indices$mean.k.[power], nrow(eigengenes),
  sum(colours == "grey"), length(colours)
))


# ---- figures -------------------------------------------------------------------------------

save_figure <- function(figure, name, width, height) {
  walk(c("png", "pdf"), \(extension) {
    ggsave(file.path(figure_dir, paste0(name, ".", extension)), figure,
      width = width, height = height, dpi = 200, bg = "white"
    )
  })
}
caption_theme <- theme(plot.caption = element_text(size = 7, colour = "grey45", hjust = 0))

save_figure(
  fit_indices |>
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
        "the 0.85 criterion and the red point the lowest power that clears it.",
        "Table: c_data/01_build_modules.xlsx, soft_threshold."
      )
    ) +
    theme_minimal(base_size = 10) +
    caption_theme,
  "01_soft_threshold", 9, 4
)

module_colours <- set_names(module_summary$module, module_summary$module)
save_figure(
  module_summary |>
    mutate(module = factor(module, levels = rev(module))) |>
    ggplot(aes(proteins, module, fill = module)) +
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
        "subject's biopsies (lme4, 1 | subject). A high value means the scored eigengene still",
        "differs between people. Table: c_data/01_build_modules.xlsx, module_summary."
      )
    ) +
    theme_minimal(base_size = 10) +
    caption_theme,
  "02_module_sizes", 7, 5
)

trajectory <- eigengene_long |>
  summarise(
    mean = mean(eigengene), se = sd(eigengene) / sqrt(n()),
    .by = c(module, arm, timepoint)
  )
save_figure(
  ggplot(trajectory, aes(timepoint, mean, colour = arm, group = arm)) +
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
    caption_theme,
  "03_trajectories", 11, 7
)

packages <- c("here", "WGCNA", "lme4", "dplyr", "purrr", "ggplot2")
versions <- tibble(
  package = packages, version = map_chr(packages, \(p) as.character(packageVersion(p)))
)
saveRDS(
  list(
    eigengenes = eigengenes, membership = membership, meta = meta,
    module_summary = module_summary, power = power, soft_threshold = fit_indices,
    provenance = list(
      created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
      inputs = manifest, packages = versions
    )
  ),
  file.path(out, "modules.rds"),
  compress = "xz"
)
writexl::write_xlsx(
  list(
    module_summary = module_summary, membership = membership,
    eigengenes = as_tibble(t(eigengenes), rownames = "sample_id"),
    soft_threshold = fit_indices, input_manifest = manifest, package_versions = versions
  ),
  file.path(out, "01_build_modules.xlsx")
)
combined <- file.path(figure_dir, "01_build_modules_figures.pdf")
pages <- setdiff(list.files(figure_dir, "[.]pdf$", full.names = TRUE), combined)
invisible(qpdf::pdf_combine(sort(pages), combined))
message("wrote modules.rds, 01_build_modules.xlsx and a ", length(pages), "-page figure PDF")

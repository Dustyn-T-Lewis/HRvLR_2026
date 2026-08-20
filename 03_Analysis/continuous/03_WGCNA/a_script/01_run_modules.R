#!/usr/bin/env Rscript
# Module construction, continuous tree. Mirrors
# categorical/03_WGCNA/a_script/setup.R's module-building block using the
# same shared engine (functions/shared_wgcna.R) on the same input matrix --
# module construction never references Group, so this build is expected to
# be numerically identical to the categorical tree's, and is run
# independently anyway so each tree stays self-contained and deletable on
# its own -- confirmed byte-identical to the categorical rebuild.
#
# Scope note: this script writes the eigengene matrix both trees' contrast
# fits and the supervised sub-stage need. The categorical tree's full
# diagnostic composite (module atlas, GO:BP ORA, module-contrast heatmap,
# baseline-eigengene LOSO prediction, rendered WGCNA.png) is not mirrored
# here yet -- deferred, not dropped.
pacman::p_load(here, dplyr, tidyr, tibble, readr)

source(here("functions", "shared_wgcna.R"))

DAT_DIR <- here("03_Analysis", "continuous", "03_WGCNA", "c_data")
dir.create(DAT_DIR, recursive = TRUE, showWarnings = FALSE)

imputed <- readRDS(here(
  "02_Normalization", "imputation", "c_data", "DAList_imputed_missforest.rds"
))
stopifnot(identical(imputed$metadata$Col_ID, colnames(imputed$data)))

abund <- as.matrix(imputed$data)
rownames(abund) <- imputed$annotation$gene

meta <- imputed$metadata |>
  transmute(sample_id = Col_ID, subject = Subject_ID)

wg <- fit_modules(abund, meta$subject)
mods <- wg$colors
me_long <- eigengene_long(wg$eigengenes, meta)

cat(sprintf(
  paste0(
    "continuous modules: power %d (R2 = %.3f, mean k = %.1f) | ",
    "%d modules, %d grey of %d\n"
  ),
  wg$power, wg$r2, wg$mean_k,
  length(setdiff(unique(mods), "grey")), sum(mods == "grey"), length(mods)
))

write_csv(
  me_long |> transmute(group_id = sub("^ME", "", module), sample_id, ME),
  file.path(DAT_DIR, "wgcna_eigengene.csv")
)

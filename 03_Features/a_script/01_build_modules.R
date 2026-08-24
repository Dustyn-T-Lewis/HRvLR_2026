# Co-expression modules on the imputed proteome, one of the three feature
# levels the association stage fits. Twelve eigengenes against fourteen
# subjects is a workable ratio where 1900 proteins is not, so this level is
# where a real effect has the best chance of clearing multiplicity.
#
# Module construction never sees a phenotype, so the modules are fixed before
# any association is tested. The engine's one deliberate deviation is
# unchanged and documented in shared_wgcna.R:
# modules are defined on within-subject-centred abundance so subject identity
# cannot drive them, then scored on raw abundance so between-subject contrasts
# survive to be tested.

pacman::p_load(here, dplyr, tibble, readr, openxlsx)

source(here("functions", "shared_wgcna.R"))

OUT_DIR <- here("03_Features", "c_data")

imputed <- readRDS(here(
  "02_Normalization", "imputation", "c_data", "DAList_imputed_missforest.rds"
))
stopifnot(identical(imputed$metadata$Col_ID, colnames(imputed$data)))

abund <- as.matrix(imputed$data)
rownames(abund) <- imputed$annotation$gene

meta <- imputed$metadata |>
  transmute(
    sample_id = Col_ID, subject = Subject_ID,
    group = factor(Group, levels = c("HR", "LR")),
    timepoint = factor(Timepoint, levels = c("T1", "T2", "T3"))
  )

set.seed(42)
wg <- fit_modules(abund, meta$subject)
me_long <- eigengene_long(wg$eigengenes, meta)

membership <- tibble::enframe(wg$colors, name = "gene", value = "module")

# How much eigengene variance subject identity still explains. The centring
# step above is what keeps this low; if it climbs, the modules have started
# encoding who the sample came from rather than what it is.
icc <- subject_variance(me_long)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

write_csv(
  me_long |> transmute(group_id = sub("^ME", "", module), sample_id, ME),
  file.path(OUT_DIR, "wgcna_eigengene.csv")
)
write.xlsx(
  list(
    module_membership = membership,
    module_eigengene = me_long,
    subject_icc = icc,
    soft_threshold = wg$sft$fitIndices
  ),
  file.path(OUT_DIR, "01_modules.xlsx")
)

message(sprintf(
  "power %d (signed R2 = %.3f, mean k = %.1f); %d modules, %d grey of %d",
  wg$power, wg$r2, wg$mean_k,
  length(setdiff(unique(wg$colors), "grey")),
  sum(wg$colors == "grey"), length(wg$colors)
))

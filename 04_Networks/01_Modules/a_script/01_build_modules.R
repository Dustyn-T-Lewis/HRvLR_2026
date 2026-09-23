# Co-expression modules: the module-by-sample eigengene matrix and each
# protein's membership, which every later network step reads.
#
# The engine is functions/shared_wgcna.R, unchanged: modules are defined on
# abundance centred within subject, so subject identity cannot drive them,
# then scored on raw abundance, so between-arm differences survive to be
# tested. Construction never sees a label or a phenotype. missForest input,
# because WGCNA needs a complete matrix.

pacman::p_load(here, dplyr, tibble, purrr, WGCNA)

source(here("functions", "shared_wgcna.R"))
source(here("functions", "association.R"))
source(here("functions", "shared_utils.R"))

OUT_DIR <- here("04_Networks", "01_Modules", "c_data")

imputed <- readRDS(here(
  "01_Preprocess", "03_Imputation", "c_data", "DAList_imputed_missforest.rds"
))
proteins <- readRDS(here(
  "02_Proteins", "01_Differential", "c_data", "proteins.rds"
))
abund <- as.matrix(imputed$data)
meta <- proteins$meta
stopifnot(
  identical(colnames(abund), meta$sample_id),
  identical(rownames(abund), proteins$annotation$uniprot_id)
)

set.seed(42)
wg <- fit_modules(abund, meta$subject)

eigengenes <- t(as.matrix(wg$eigengenes))
rownames(eigengenes) <- sub("^ME", "", rownames(eigengenes))
genes <- setNames(proteins$annotation$gene, proteins$annotation$uniprot_id)

kme <- WGCNA::signedKME(
  t(abund), wg$eigengenes,
  corFnc = "bicor",
  corOptions = "maxPOutliers = 0.05, pearsonFallback = 'individual'"
)
colnames(kme) <- sub("^kME", "", colnames(kme))
membership <- tibble(
  uniprot_id = names(wg$colors), gene = genes[names(wg$colors)],
  module = unname(wg$colors)
) |>
  mutate(kme = map2_dbl(
    .data$uniprot_id, .data$module,
    \(id, m) if (m == "grey") NA_real_ else kme[id, m]
  ))

icc <- subject_variance(eigengene_long(wg$eigengenes, meta)) |>
  mutate(module = sub("^ME", "", .data$module))

modules <- list(
  eigengenes = eigengenes,
  membership = membership,
  meta = meta,
  icc = icc,
  power = wg$power, r2 = wg$r2, mean_k = wg$mean_k,
  soft_threshold = wg$sft$fitIndices,
  windows = set_names(c("T1", "training", "acute")) |>
    map(\(w) subject_window(eigengenes, w, meta))
)

clear_dir(OUT_DIR)
saveRDS(modules, file.path(OUT_DIR, "modules.rds"))
openxlsx::write.xlsx(list(
  membership = membership,
  eigengenes = as_tibble(t(eigengenes), rownames = "sample_id"),
  subject_icc = icc,
  soft_threshold = wg$sft$fitIndices
), file.path(OUT_DIR, "modules.xlsx"))

message(sprintf(
  "power %d (signed R2 = %.3f, mean k = %.1f); %d modules, %d grey of %d",
  wg$power, wg$r2, wg$mean_k, nrow(eigengenes),
  sum(wg$colors == "grey"), length(wg$colors)
))

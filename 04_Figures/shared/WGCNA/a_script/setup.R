# F04 setup: read the imputed matrix, the stage-00 phenotype table, and the limma fit from
# 03_DEP, then build the modules. Panels and 01_run_modules.R source this first. Writes nothing.
pacman::p_load(here, dplyr, tidyr, tibble, readr, purrr)

source(here("functions", "shared_style.R"))
source(here("functions", "shared_pathway_utils.R"))
source(here("functions", "shared_wgcna.R"))
source(here("04_Figures", "F04_association", "WGCNA", "a_script", "enrichment.R"))
source(here("03_DEP", "contrasts.R"))

RPT_DIR <- here("04_Figures", "F04_association", "WGCNA", "b_reports")
DAT_DIR <- here("04_Figures", "F04_association", "WGCNA", "c_data")

if (!exists("F04_PANELS")) F04_PANELS <- list()
if (!exists("F04_AUDIT")) F04_AUDIT <- list()

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

pheno <- read_csv(here("00_input", "c_data", "phenotype.csv"), show_col_types = FALSE)
ADAPTATION_TRAITS <- setdiff(names(pheno), c("subject", "group_arm"))

TRAIT_LABELS <- c(
  comp_hypertrophy = "Composite hypertrophy",
  d_fcsa_I = "Δ fCSA I", d_fcsa_II = "Δ fCSA II", d_mcsa = "Δ mCSA",
  d_1rm_legpress = "Δ 1RM leg-press", d_1rm_ext = "Δ 1RM leg-ext"
)
stopifnot(setequal(names(TRAIT_LABELS), ADAPTATION_TRAITS))

# The DEP fit, read not refitted: 03_DEP already estimated the design, the contrasts, and
# the within-subject correlation. fry reuses all three rather than re-deriving them.
dep_fit <- readRDS(here("03_DEP", "a_non_imputed", "c_data", "01_limma_DAList.rds"))
DESIGN <- dep_fit$design$design_matrix
CONTRAST_MATRIX <- dep_fit$design$contrast_matrix
CORRELATION <- dep_fit$eBayes_fit$correlation

# The block is the random factor 03_DEP actually fitted (subject), named on its own metadata.
BLOCK <- dep_fit$metadata[[dep_fit$design$random_factor]]

# fry takes the design rows positionally, so the abundance columns must be in the design's
# order. They happen to agree today; assert it rather than trust it, because a silent
# mismatch would pair every sample with another sample's design row.
DEP_SAMPLES <- rownames(DESIGN)
stopifnot(
  setequal(DEP_SAMPLES, colnames(abund)),
  length(BLOCK) == nrow(DESIGN)
)
abund_dep <- abund[, DEP_SAMPLES, drop = FALSE]
stopifnot(identical(colnames(abund_dep), DEP_SAMPLES))

wg <- fit_modules(abund, meta$subject)
mods <- wg$colors
me_long <- eigengene_long(wg$eigengenes, meta)

cat(sprintf(
  "F04 modules: power %d (signed R2 = %.3f, mean k = %.1f) | %d modules, %d grey of %d\n",
  wg$power, wg$r2, wg$mean_k,
  length(setdiff(unique(mods), "grey")), sum(mods == "grey"), length(mods)
))

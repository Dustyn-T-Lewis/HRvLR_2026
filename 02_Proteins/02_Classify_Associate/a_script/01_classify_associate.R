# Protein-level classification and phenotype association.
#
# Reads the abundance matrix without imputation, the same one the contrasts
# were fitted on. A protein missing in a subject drops out of that subject's
# pair rather than being filled in, so a small-n AUC is honest about its n.

pacman::p_load(here, dplyr)

source(here("functions", "classify.R"))
source(here("functions", "shared_utils.R"))

OUT_DIR <- here("02_Proteins", "02_Classify_Associate", "c_data")

proteins <- readRDS(here(
  "02_Proteins", "01_Differential", "c_data", "proteins.rds"
))
genes <- setNames(proteins$annotation$gene, proteins$annotation$uniprot_id)

screens <- screen_level(proteins$abund, proteins$meta)
screens[c("classify", "associate")] <- lapply(
  screens[c("classify", "associate")],
  \(d) mutate(d, gene = genes[.data$feature], .after = "feature")
)

clear_dir(OUT_DIR)
saveRDS(screens, file.path(OUT_DIR, "protein_screens.rds"))
write_workbook(file.path(OUT_DIR, "protein_screens.xlsx"), screens)

print(as.data.frame(screens$chance_classify), row.names = FALSE, digits = 3)

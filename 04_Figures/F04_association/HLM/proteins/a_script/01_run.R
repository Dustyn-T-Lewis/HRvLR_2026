# Per-protein association with each continuous training response, trajectory
# view: all three timepoints kept, random intercept per subject, moderated
# with a mixed model interaction F-test. Descriptive, n = 16, exploratory.
pacman::p_load(here, dplyr, readr, purrr, parallel, openxlsx)

source(here("functions", "shared_style.R"))
source(here("04_Figures", "functions", "f04_association.R"))

leaf <- here("04_Figures", "F04_association", "HLM", "proteins")
data_dir <- file.path(leaf, "c_data")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

dal <- readRDS(here("02_Normalization", "c_data", "DAList_normalized.rds"))
meta <- dal$metadata
pheno <- read_csv(here("00_input", "c_data", "phenotype.csv"), show_col_types = FALSE)
pheno$subject <- meta$Subject_ID[match(pheno$subject, meta$Subject_ID)]
traits <- intersect(
  c(
    "comp_hypertrophy", "d_fcsa_I", "d_fcsa_II",
    "d_mcsa", "d_1rm_legpress", "d_1rm_ext"
  ),
  names(pheno)
)

cores <- as.integer(Sys.getenv(
  "ASSOC_CORES", as.character(max(1L, parallel::detectCores() - 2L))
))

gene_of <- setNames(dal$annotation$gene, rownames(dal$data))
protein_assoc <- associate_hlm(dal$data, meta, pheno, traits, cores = cores) |>
  mutate(gene = gene_of[feature], .after = feature)

wb <- createWorkbook()
addWorksheet(wb, "protein_association_hlm")
writeData(wb, "protein_association_hlm", protein_assoc)
saveWorkbook(wb, file.path(data_dir, "proteins_source_data.xlsx"), overwrite = TRUE)

message(sprintf(
  "protein HLM: %d finite interaction p, %d NA (sparse guard), %d BH < 0.05",
  sum(!is.na(protein_assoc$p)), sum(is.na(protein_assoc$p)),
  sum(protein_assoc$bh < 0.05, na.rm = TRUE)
))

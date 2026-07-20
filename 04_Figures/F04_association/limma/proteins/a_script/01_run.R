# Per-protein association with each continuous training response, snapshot view:
# each subject collapsed to one value per phase (baseline level, training and
# acute change), moderated with limma. Descriptive, n = 16, exploratory.
pacman::p_load(here, dplyr, readr, purrr, ggplot2, openxlsx)


source(here("functions", "shared_style.R"))
source(here("04_Figures", "functions", "f04_association.R"))
source(here("04_Figures", "functions", "f04_plots.R"))

leaf <- here("04_Figures", "F04_association", "limma", "proteins")
data_dir <- file.path(leaf, "c_data")
report_dir <- file.path(leaf, "b_reports")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

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

gene_of <- setNames(dal$annotation$gene, rownames(dal$data))
protein_assoc <- associate_limma(dal$data, meta, pheno, traits) |>
  mutate(gene = gene_of[feature], .after = feature)
wb <- createWorkbook()
addWorksheet(wb, "protein_association")
writeData(wb, "protein_association", protein_assoc)
saveWorkbook(wb, file.path(data_dir, "proteins_source_data.xlsx"), overwrite = TRUE)

lead_trait <- if ("d_mcsa" %in% traits) "d_mcsa" else traits[1]
volc <- assoc_volcano(
  dplyr::filter(protein_assoc, trait == lead_trait, phase == "baseline"),
  trait_label = paste("baseline proteome vs", lead_trait)
)
ggsave(
  file.path(report_dir, "protein_volcano_mcsa.pdf"), volc,
  width = 6, height = 5
)

message(sprintf(
  "protein association: %d BH < 0.05 of %d rows",
  sum(protein_assoc$bh < 0.05, na.rm = TRUE), nrow(protein_assoc)
))

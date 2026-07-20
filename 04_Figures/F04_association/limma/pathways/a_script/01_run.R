# Per-pathway (singscore) association with each continuous training response,
# same snapshot design as the proteins leaf. Proteins collapse to one row per
# gene (highest mean), keep genes detected in every sample so the rank-based
# score is well defined, then score against the shared pathway collection.
pacman::p_load(here, dplyr, readr, tibble, purrr, ggplot2, openxlsx)

set.seed(42)

source(here("functions", "shared_style.R"))
source(here("functions", "shared_pathway_utils.R"))
source(here("functions", "shared_singscore.R"))
source(here("04_Figures", "functions", "f04_association.R"))
source(here("04_Figures", "functions", "f04_plots.R"))

leaf <- here("04_Figures", "F04_association", "limma", "pathways")
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
gene_mat <- as.data.frame(dal$data) |>
  mutate(gene = gene_of[rownames(dal$data)]) |>
  filter(!is.na(gene)) |>
  mutate(.m = rowMeans(across(where(is.numeric)), na.rm = TRUE)) |>
  group_by(gene) |>
  slice_max(.m, n = 1, with_ties = FALSE) |>
  ungroup()
mat <- as.matrix(gene_mat[, meta$Col_ID])
rownames(mat) <- gene_mat$gene
mat <- mat[rowSums(is.na(mat)) == 0, , drop = FALSE]
message(sprintf(
  "singscore matrix: %d complete genes x %d samples", nrow(mat), ncol(mat)
))

gene_sets <- build_pathway_collection()
scores <- score_singscore(mat, gene_sets, min_size = 5)
message(sprintf("scored %d pathways", nrow(scores)))

pathway_scores <- scores |>
  as.data.frame() |>
  rownames_to_column("pathway")

pathway_assoc <- associate_limma(scores, meta, pheno, traits)

wb <- createWorkbook()
addWorksheet(wb, "pathway_association")
writeData(wb, "pathway_association", pathway_assoc)
addWorksheet(wb, "pathway_scores")
writeData(wb, "pathway_scores", pathway_scores)
saveWorkbook(wb, file.path(data_dir, "pathways_source_data.xlsx"), overwrite = TRUE)

lead_trait <- if ("d_mcsa" %in% traits) "d_mcsa" else traits[1]
volc <- assoc_volcano(
  dplyr::filter(pathway_assoc, trait == lead_trait, phase == "baseline"),
  trait_label = paste("baseline pathway scores vs", lead_trait)
)
ggsave(
  file.path(report_dir, "pathway_volcano_mcsa.pdf"), volc,
  width = 6, height = 5
)

message(sprintf(
  "pathway association: %d BH < 0.05 of %d rows",
  sum(pathway_assoc$bh < 0.05, na.rm = TRUE), nrow(pathway_assoc)
))

# Per-pathway (singscore) association with each continuous training response,
# same trajectory design as the proteins leaf. Proteins collapse to one row
# per gene (highest mean), keep genes detected in every sample so the
# rank-based score is well defined, then score against the shared pathway
# collection.
pacman::p_load(here, dplyr, readr, tibble, purrr, parallel, openxlsx)

source(here("functions", "shared_style.R"))
source(here("functions", "shared_pathway_utils.R"))
source(here("functions", "shared_singscore.R"))
source(here("04_Figures", "functions", "f04_association.R"))

leaf <- here("04_Figures", "F04_association", "HLM", "pathways")
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

pathway_assoc <- associate_hlm(scores, meta, pheno, traits, cores = cores)

wb <- createWorkbook()
addWorksheet(wb, "pathway_association_hlm")
writeData(wb, "pathway_association_hlm", pathway_assoc)
addWorksheet(wb, "pathway_scores")
writeData(wb, "pathway_scores", pathway_scores)
saveWorkbook(wb, file.path(data_dir, "pathways_source_data.xlsx"), overwrite = TRUE)

message(sprintf(
  "pathway HLM: %d finite interaction p, %d NA (sparse guard), %d BH < 0.05",
  sum(!is.na(pathway_assoc$p)), sum(is.na(pathway_assoc$p)),
  sum(pathway_assoc$bh < 0.05, na.rm = TRUE)
))

# Phenotype mapping: which proteins and pathways associate with the continuous
# training responses (mCSA, strength, fibre CSA). Association only, per the
# small-n methods review in docs/phenotype_mapping_methods_review.md; prediction
# is out of scope here and belongs in a clearly-null supplement. Writes tidy
# association tables and the singscore matrix to c_data, and lead figures to
# b_reports.

pacman::p_load(here, dplyr, readr, tibble, purrr, ggplot2)

unit <- here("04_Figures", "extras", "association")
source(here("04_Figures", "functions", "style.R"))
source(here("04_Figures", "functions", "pathway_utils.R"))
source(file.path(unit, "a_script", "functions", "associate.R"))
source(file.path(unit, "a_script", "functions", "plots.R"))

data_dir <- file.path(unit, "c_data")
report_dir <- file.path(unit, "b_reports")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

dal <- readRDS(here("02_Normalization", "c_data", "DAList_normalized.rds"))
meta <- dal$metadata
pheno <- readr::read_csv(
  here("04_Figures", "F04", "c_data", "phenotype.csv"),
  show_col_types = FALSE
)
traits <- intersect(
  c(
    "comp_hypertrophy", "d_fcsa_I", "d_fcsa_II",
    "d_mcsa", "d_1rm_legpress", "d_1rm_ext"
  ),
  names(pheno)
)

# Align phenotype subjects to the DAList subject labels.
pheno$subject <- meta$Subject_ID[match(pheno$subject, meta$Subject_ID)]
if (all(is.na(pheno$subject))) {
  stop("phenotype `subject` does not match DAList Subject_ID; check the key")
}

phases <- c("baseline", "training", "acute")

# Protein-level association on the NA-aware normalized matrix.
gene_of <- setNames(dal$annotation$gene, rownames(dal$data))
protein_assoc <- purrr::map_dfr(
  phases, \(ph) associate_traits(dal$data, meta, pheno, traits, ph)
) |>
  dplyr::mutate(gene = gene_of[feature], .after = feature)
readr::write_csv(protein_assoc, file.path(data_dir, "protein_association.csv"))

# Pathway scoring: collapse to one row per gene (highest mean), keep proteins
# detected in every sample so the rank-based score is well defined.
gene_mat <- as.data.frame(dal$data) |>
  dplyr::mutate(gene = gene_of[rownames(dal$data)]) |>
  dplyr::filter(!is.na(gene)) |>
  dplyr::mutate(
    .m = rowMeans(dplyr::across(dplyr::where(is.numeric)), na.rm = TRUE)
  ) |>
  dplyr::group_by(gene) |>
  dplyr::slice_max(.m, n = 1, with_ties = FALSE) |>
  dplyr::ungroup()
mat <- as.matrix(gene_mat[, meta$Col_ID])
rownames(mat) <- gene_mat$gene
mat <- mat[rowSums(is.na(mat)) == 0, , drop = FALSE]
message(sprintf(
  "singscore matrix: %d complete genes x %d samples", nrow(mat), ncol(mat)
))

gene_sets <- build_pathway_collection()
scores <- score_pathways(mat, gene_sets)
message(sprintf("scored %d pathways", nrow(scores)))
scores |>
  as.data.frame() |>
  tibble::rownames_to_column("pathway") |>
  readr::write_csv(file.path(data_dir, "pathway_scores.csv"))

pathway_assoc <- purrr::map_dfr(
  phases, \(ph) associate_traits(scores, meta, pheno, traits, ph)
)
readr::write_csv(pathway_assoc, file.path(data_dir, "pathway_association.csv"))

# Lead figures: protein volcano for mCSA, and a singscore rank-density for the
# top mCSA-associated pathway.
lead_trait <- if ("d_mcsa" %in% traits) "d_mcsa" else traits[1]
volc <- assoc_volcano(
  dplyr::filter(protein_assoc, trait == lead_trait, phase == "baseline"),
  trait_label = paste("baseline proteome vs", lead_trait)
)
ggsave(
  file.path(report_dir, "protein_volcano_mcsa.pdf"), volc,
  width = 6, height = 5
)

top_path <- pathway_assoc |>
  dplyr::filter(trait == lead_trait) |>
  dplyr::slice_min(bh, n = 1, with_ties = FALSE) |>
  dplyr::pull(feature)
rank_fig <- singscore_rank_density(mat, gene_sets[[top_path]], meta$Col_ID[1])
ggsave(
  file.path(report_dir, "singscore_rank_density.pdf"), rank_fig,
  width = 6, height = 4
)

message("phenotype_mapping done -> ", data_dir)

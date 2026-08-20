# The three feature matrices the nine contrasts are fitted on, all keyed by the
# 45 sample ids. Proteins are the non-imputed cycloess matrix stage 03 fits;
# modules and pathways are derived from the missForest arm, an asymmetry that
# every figure built on them has to state.

pacman::p_load(here, dplyr, tidyr, readr)

SINGSCORE_CACHE <- c(
  "03_Analysis", "categorical", "02_Pathways", "c_data", "singscore_scores.rds"
)

LEVEL_DIR <- c(
  pathways = "02_Pathways", modules = "03_WGCNA", proteins = "01_Proteins"
)

protein_matrix <- function() {
  df <- read_csv(here("02_Normalization", "c_data", "normalized.csv"),
    show_col_types = FALSE
  )
  ann_cols <- c("uniprot_id", "protein", "gene", "description")
  mat <- as.matrix(df[, setdiff(names(df), ann_cols)])
  rownames(mat) <- df$uniprot_id
  mat
}

module_matrix <- function() {
  eigen_csv <- here(
    "03_Analysis", "categorical", "03_WGCNA", "c_data", "wgcna_eigengene.csv"
  )
  wide <- read_csv(eigen_csv, show_col_types = FALSE) |>
    pivot_wider(names_from = "sample_id", values_from = "ME") |>
    as.data.frame()
  mat <- as.matrix(wide[, setdiff(names(wide), "group_id")])
  rownames(mat) <- paste0("ME_", wide$group_id)
  mat
}

pathway_matrix <- function() {
  readRDS(do.call(here, as.list(SINGSCORE_CACHE)))
}

# Sets are filtered to SET_FLOOR detected members before scoring, matching what
# fgsea and fry apply. The counts still travel with the row so a reader can see
# how much of each set was measured.
pathway_coverage <- function(gene_sets, detected_genes) {
  tibble(
    feature = names(gene_sets),
    n_annotated = lengths(gene_sets),
    n_detected = vapply(
      gene_sets, function(g) sum(g %in% detected_genes), integer(1)
    )
  ) |>
    mutate(coverage = .data$n_detected / .data$n_annotated)
}

detected_genes <- function() {
  df <- read_csv(here("02_Normalization", "c_data", "normalized.csv"),
    show_col_types = FALSE
  )
  unique(df$gene[!is.na(df$gene) & df$gene != ""])
}

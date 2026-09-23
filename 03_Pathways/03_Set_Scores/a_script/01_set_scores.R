# Per-sample pathway scores: the set-by-sample matrix the pathway screens and
# the packet read.
#
# singscore ranks within each sample, so a score depends only on that sample's
# own protein ranks and carries nothing from the rest of the cohort or from any
# label. Scored on the missForest matrix because a rank needs every protein
# present; the sets are the same 15-to-500 detected-member sets fry tested.

pacman::p_load(here, dplyr, purrr, singscore, GSEABase)

source(here("functions", "association.R"))

OUT_DIR <- here("03_Pathways", "03_Set_Scores", "c_data")

gs <- readRDS(here("03_Pathways", "01_Gene_Sets", "c_data", "gene_sets.rds"))
proteins <- readRDS(here(
  "02_Proteins", "01_Differential", "c_data", "proteins.rds"
))
imputed <- readRDS(here(
  "01_Preprocess", "03_Imputation", "c_data", "DAList_imputed_missforest.rds"
))
abund <- as.matrix(imputed$data)
stopifnot(identical(colnames(abund), proteins$meta$sample_id))

collection <- GeneSetCollection(
  imap(gs$sets, \(ids, name) GeneSet(ids, setName = name))
)
scores <- multiScore(rankGenes(abund), upSetColc = collection)$Scores
stopifnot(identical(rownames(scores), gs$catalog$set))

set_scores <- list(
  scores = scores,
  meta = proteins$meta,
  catalog = gs$catalog,
  windows = set_names(c("T1", "training", "acute")) |>
    map(\(w) subject_window(scores, w, proteins$meta))
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
saveRDS(set_scores, file.path(OUT_DIR, "set_scores.rds"))
readr::write_csv(
  tibble::as_tibble(scores, rownames = "set"),
  file.path(OUT_DIR, "set_scores.csv")
)
message(sprintf("%d sets x %d samples", nrow(scores), ncol(scores)))

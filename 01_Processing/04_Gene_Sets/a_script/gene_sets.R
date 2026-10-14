# Read one MSigDB release, map the protein matrix onto gene symbols, and keep the sets large enough
# to test. Writes the set list every set-level step in 01_Processing and the analysis folders read,
# so membership is decided once.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tibble)
  library(purrr)
})
source(here("00_Setup", "Shared_functions", "analysis_reporting_helpers.R"))

out <- here("01_Processing", "04_Gene_Sets", "c_data")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

inputs <- c(
  proteins = "01_Processing/02_Normalisation/c_data/DAList_normalized.rds",
  msigdb = "01_Processing/04_Gene_Sets/references/msigdb_2026.1.Hs.rds",
  goslim = "01_Processing/04_Gene_Sets/references/goslim_generic.obo"
)
paths <- map_chr(inputs, here)
if (!all(file.exists(paths))) {
  stop(
    "Run 01_Processing/02_Normalisation first. Missing: ",
    paste(inputs[!file.exists(paths)], collapse = ", ")
  )
}
proteins <- readRDS(paths[["proteins"]])

# Membership shifts between MSigDB releases, so one release (2026.1.Hs) is read from a saved
# msigdbr snapshot rather than from the installed package.
msigdb <- readRDS(paths[["msigdb"]])
collections <- msigdb$collections
membership <- rename(msigdb$membership, collection = database)

# Set tests need one row per gene. A protein with no symbol or several is dropped rather than
# split. Where proteins share a symbol, the one observed in the most samples represents it.
# Decided once, reading no fold change.
protein_map <- proteins$annotation |>
  as_tibble() |>
  transmute(
    uniprot_id, n_seq,
    n_obs = rowSums(!is.na(proteins$data))[uniprot_id],
    gene = trimws(gene),
    mapping_status = case_when(
      is.na(gene) | !nzchar(gene) ~ "missing_symbol",
      grepl("[;,|]", gene) ~ "multiple_symbols",
      TRUE ~ "candidate"
    ),
    gene = if_else(mapping_status == "candidate", gene, NA_character_)
  )
representatives <- protein_map |>
  filter(mapping_status == "candidate") |>
  arrange(gene, desc(n_obs), desc(n_seq), uniprot_id) |>
  distinct(gene, .keep_all = TRUE) |>
  pull(uniprot_id)
protein_map <- protein_map |>
  mutate(
    selected = uniprot_id %in% representatives,
    mapping_status = case_when(
      selected ~ "representative",
      mapping_status == "candidate" ~ "duplicate_gene",
      TRUE ~ mapping_status
    ),
    # Accessions distinguish duplicate symbols in protein-level plot labels.
    label = if_else(
      is.na(gene), uniprot_id,
      if_else(duplicated(gene) | duplicated(gene, fromLast = TRUE),
        paste0(gene, " (", uniprot_id, ")"), gene
      )
    )
  )
gene_map <- filter(protein_map, selected)
gene_universe <- gene_map$gene
stopifnot(!anyDuplicated(protein_map$label))
mapping_summary <- count(protein_map, mapping_status, name = "n_proteins")
message("measured gene universe: ", length(gene_universe))

# Size is the only pre-test filter, and the measured bar counts what this experiment detected.
# Nothing is dropped for overlapping another set or for its name.
set_members <- distinct(membership, set_id, gene)
sets_full <- split(set_members$gene, set_members$set_id)
sets_measured <- map(sets_full, intersect, y = gene_universe)

set_catalog <- membership |>
  distinct(set_id, collection, pathway, source_id, description) |>
  mutate(
    source_size = lengths(sets_full)[set_id],
    measured_size = lengths(sets_measured)[set_id],
    qualifies = source_size >= 15 & source_size <= 500 & measured_size >= 15
  )
stopifnot(!anyDuplicated(set_catalog$set_id))
sets <- sets_measured[set_catalog$set_id[set_catalog$qualifies]]
if (!length(sets)) stop("No gene sets passed the size filters.")

# GO Slim sets come from the GO Consortium's generic slim.
slim_offspring <- AnnotationDbi::mget(
  GSEABase::ids(GSEABase::getOBOCollection(paths[["goslim"]])), GO.db::GOBPOFFSPRING,
  ifnotfound = NA
)
slim_offspring <- slim_offspring[!is.na(slim_offspring)]

# Each slim term is also a set: every measured gene annotated to it or to any term beneath it,
# taken from the full MSigDB membership so a gene is not lost when its GO:BP set falls outside
# the size filter. The 15-to-500 rule reads measured size, since a slim term is broad by design.
go_genes <- membership |>
  filter(collection == "GOBP") |>
  with(split(gene, source_id))
slim_sets <- imap(slim_offspring, function(descendants, slim_id) {
  covered <- intersect(c(slim_id, descendants), names(go_genes))
  sort(intersect(unique(unlist(go_genes[covered], use.names = FALSE)), gene_universe))
})
slim_catalog <- tibble(
  theme_id = names(slim_sets),
  pathway = unname(AnnotationDbi::Term(GO.db::GOTERM[theme_id])),
  measured_size = lengths(slim_sets)
) |>
  transmute(
    set_id = paste("GO_Slim", toupper(gsub("[^A-Za-z0-9]+", "_", pathway)), sep = "|"),
    collection = "GO_Slim", pathway, source_id = theme_id,
    description = "GO Slim term: every measured gene under it in the GO:BP hierarchy",
    source_size = measured_size, measured_size,
    qualifies = measured_size >= 15 & measured_size <= 500
  )
names(slim_sets) <- slim_catalog$set_id
set_catalog <- bind_rows(set_catalog, slim_catalog)
sets <- c(sets, slim_sets[slim_catalog$set_id[slim_catalog$qualifies]])
stopifnot(!anyDuplicated(set_catalog$set_id), !anyDuplicated(names(sets)))
message("GO Slim sets: ", sum(slim_catalog$qualifies), " of ", nrow(slim_catalog), " testable")

collection_summary <- set_catalog |>
  summarise(
    n_in_source = n(), n_qualifying = sum(qualifies),
    median_measured = median(measured_size[qualifies]), .by = collection
  ) |>
  arrange(match(collection, c(collections, "GO_Slim")))
message("qualifying sets: ", length(sets))

saveRDS(
  list(
    sets = sets, set_catalog = set_catalog, protein_map = protein_map,
    gene_universe = gene_universe
  ),
  file.path(out, "gene_sets.rds"),
  compress = "xz"
)
sheets <- list(
  collection_summary = collection_summary,
  set_catalog = set_catalog,
  protein_gene_map = protein_map,
  mapping_summary = mapping_summary
)
holds <- c(
  "Sets per collection in the source, sets that qualify, median measured size.",
  "Every set with its source and measured size and whether it qualifies (15 to 500).",
  "Every protein, its gene symbol, and the representative chosen per symbol.",
  "Proteins per mapping status."
)
write_workbook(sheets, holds, file.path(out, "gene_sets.xlsx"), inputs, map_chr(inputs, here))

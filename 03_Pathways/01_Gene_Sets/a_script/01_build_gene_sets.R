# The gene-set universe every later pathway and module test reads.
#
# Three MSigDB collections at a pinned release, plus one set per GO-Slim term.
# Sets are rewritten from gene symbols to the protein ids the fit uses, each
# symbol to one protein: where two proteins share a symbol, the one observed
# in more samples carries it. A set is kept when 15 to 500 of its members were
# detected, the same floor every downstream test applies.
#
# Every GO:BP set also gets a theme, its most specific GO-Slim ancestor, which
# is how pathways are grouped in the packet. Hallmark and Reactome sets carry
# no GO identifier and so no theme.

pacman::p_load(here, dplyr, tibble, purrr, msigdbr)

source(here("functions", "shared_pathway_utils.R"))
source(here("functions", "shared_utils.R"))

MSIGDB_RELEASE <- "2026.1.Hs"
OUT_DIR <- here("03_Pathways", "01_Gene_Sets", "c_data")

proteins <- readRDS(here(
  "02_Proteins", "01_Differential", "c_data", "proteins.rds"
))
protein_map <- proteins$annotation |>
  mutate(n_obs = rowSums(!is.na(proteins$abund))) |>
  filter(!is.na(.data$gene), .data$gene != "") |>
  slice_max(.data$n_obs, n = 1, by = "gene", with_ties = FALSE) |>
  select("gene", "uniprot_id")

COLLECTIONS <- tribble(
  ~collection, ~msig_collection, ~subcollection,
  "Hallmark",  "H",              NA,
  "Reactome",  "C2",             "CP:REACTOME",
  "GO:BP",     "C5",             "GO:BP"
)

msig <- pmap(COLLECTIONS, function(collection, msig_collection, subcollection) {
  args <- list(species = "Homo sapiens", collection = msig_collection)
  if (!is.na(subcollection)) args$subcollection <- subcollection
  do.call(msigdbr::msigdbr, args) |>
    transmute(
      collection = collection, set = .data$gs_name,
      source_id = .data$gs_exact_source, gene = .data$gene_symbol,
      db_version = .data$db_version
    )
}) |>
  list_rbind()
stopifnot(
  "MSigDB release moved; update MSIGDB_RELEASE deliberately" =
    all(msig$db_version == MSIGDB_RELEASE)
)

slim <- read_goslim_bp()
slim_long <- goslim_sets(slim) |>
  enframe("set", "gene") |>
  tidyr::unnest("gene") |>
  mutate(collection = "GO Slim", source_id = NA_character_)

long <- bind_rows(select(msig, -"db_version"), slim_long) |>
  distinct(.data$collection, .data$set, .data$gene, .keep_all = TRUE)

catalog <- long |>
  summarise(
    source_id = first(.data$source_id),
    size_source = n(),
    size_detected = sum(.data$gene %in% protein_map$gene),
    .by = c("collection", "set")
  ) |>
  filter(between(.data$size_detected, SET_FLOOR, SET_CEILING)) |>
  mutate(theme = if_else(
    .data$collection == "GO:BP", goslim_theme(.data$source_id, slim), NA
  ))

sets <- long |>
  semi_join(catalog, by = c("collection", "set")) |>
  inner_join(protein_map, by = "gene") |>
  (\(d) split(d$uniprot_id, d$set))()
stopifnot(setequal(names(sets), catalog$set))

gene_sets <- list(
  sets = sets[catalog$set],
  catalog = catalog,
  protein_map = protein_map,
  universe = proteins$annotation$uniprot_id,
  msigdb_release = MSIGDB_RELEASE,
  goslim_release = sub(".*releases/([0-9-]+)/.*", "\\1", grep(
    "^data-version", readLines(here("00_input", "goslim_generic.obo")),
    value = TRUE
  ))
)

clear_dir(OUT_DIR)
saveRDS(gene_sets, file.path(OUT_DIR, "gene_sets.rds"))
write_workbook(file.path(OUT_DIR, "gene_sets.xlsx"), list(
  catalog = catalog,
  protein_map = protein_map,
  themes = count(filter(catalog, !is.na(.data$theme)), .data$theme, sort = TRUE)
))

print(count(catalog, .data$collection))
message(sprintf(
  "GO:BP sets with a theme: %d of %d",
  sum(!is.na(catalog$theme)), sum(catalog$collection == "GO:BP")
))

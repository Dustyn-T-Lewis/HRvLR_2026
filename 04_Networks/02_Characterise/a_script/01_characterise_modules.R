# What each module is: its enriched gene sets, its hub proteins, its GO-Slim
# theme, and whether its members interact more than chance in STRING.
#
# Enrichment is clusterProfiler's hypergeometric ORA against the pathway
# stage's own gene sets, with the 1900 detected proteins as universe; testing
# against the genome would reward any module for being muscle. BH runs within
# each module.
#
# The STRING check is STRINGdb's PPI enrichment at combined score >= 700, with
# the background set to the detected proteins STRING can map, so a module is
# compared with this proteome rather than the genome. STRINGdb's expected
# edge count accounts for each member's degree, which a plain label shuffle
# does not. Read from the local v12 files in 00_input/downloads.

pacman::p_load(here, dplyr, tidyr, tibble, purrr, clusterProfiler, STRINGdb)

source(here("functions", "shared_pathway_utils.R"))
source(here("functions", "shared_utils.R"))

OUT_DIR <- here("04_Networks", "02_Characterise", "c_data")
STRING_MIN_SCORE <- 700L
N_HUBS <- 10L

modules <- readRDS(here("04_Networks", "01_Modules", "c_data", "modules.rds"))
gs <- readRDS(here("03_Pathways", "01_Gene_Sets", "c_data", "gene_sets.rds"))
members <- filter(modules$membership, .data$module != "grey")

term2gene <- enframe(gs$sets, "set", "uniprot_id") |> unnest("uniprot_id")
ora_fit <- compareCluster(
  split(members$uniprot_id, members$module),
  fun = "enricher", TERM2GENE = term2gene, universe = gs$universe,
  minGSSize = SET_FLOOR, maxGSSize = SET_CEILING,
  pvalueCutoff = 1, qvalueCutoff = 1
)
ora <- as_tibble(ora_fit) |>
  transmute(
    module = as.character(.data$Cluster), set = .data$ID,
    gene_ratio = .data$GeneRatio, bg_ratio = .data$BgRatio,
    count = .data$Count, p = .data$pvalue, bh = .data$p.adjust
  ) |>
  left_join(select(gs$catalog, "set", "collection", "theme"), by = "set")

hubs <- members |>
  slice_max(.data$kme, n = N_HUBS, by = "module", with_ties = FALSE)

module_theme <- ora |>
  filter(.data$bh < 0.05, !is.na(.data$theme)) |>
  count(.data$module, .data$theme, name = "n_sets") |>
  slice_max(.data$n_sets, n = 1, by = "module", with_ties = FALSE)

string_db <- STRINGdb$new(
  version = "12.0", species = 9606, score_threshold = STRING_MIN_SCORE,
  network_type = "full", input_directory = here("00_input", "downloads")
)
mapped <- string_db$map(
  as.data.frame(modules$membership), "uniprot_id",
  removeUnmappedRows = TRUE
)
string_db$set_background(mapped$STRING_id)
string_check <- mapped |>
  filter(.data$module != "grey") |>
  summarise(n_mapped = n(), ids = list(.data$STRING_id), .by = "module") |>
  mutate(
    enrichment = map(.data$ids, string_db$get_ppi_enrichment),
    edges = map_dbl(.data$enrichment, "edges"),
    expected = map_dbl(.data$enrichment, "lambda"),
    ratio = .data$edges / .data$expected,
    p = map_dbl(.data$enrichment, "enrichment"),
    bh = p.adjust(.data$p, "BH")
  ) |>
  select(-"ids", -"enrichment")

top_set <- ora |>
  slice_min(.data$p, n = 1, by = "module", with_ties = FALSE) |>
  select("module", top_set = "set", top_set_bh = "bh")
labels <- count(members, .data$module, name = "n_proteins") |>
  left_join(top_set, by = "module") |>
  left_join(select(module_theme, "module", "theme"), by = "module") |>
  left_join(
    summarise(hubs, hubs = paste(.data$gene, collapse = ", "), .by = "module"),
    by = "module"
  ) |>
  left_join(
    select(string_check, "module", string_ratio = "ratio", string_bh = "bh"),
    by = "module"
  )

clear_dir(OUT_DIR)
characterisation <- list(
  labels = labels, ora = ora, hubs = hubs, string = string_check,
  ora_fit = ora_fit, string_min_score = STRING_MIN_SCORE
)
saveRDS(characterisation, file.path(OUT_DIR, "module_characterisation.rds"))
openxlsx::write.xlsx(
  characterisation[c("labels", "ora", "hubs", "string")],
  file.path(OUT_DIR, "module_characterisation.xlsx")
)

print(as.data.frame(select(string_check, -"n_mapped")), digits = 3)

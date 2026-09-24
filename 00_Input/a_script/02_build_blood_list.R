# Collapse the five published red-cell proteome sources in RBC_proteins.xlsx into one lookup:
# one row per protein, the sources listing it, and how many agree. Mature erythrocytes have
# no nucleus, so no transcript atlas covers their proteome; these mass-spectrometry references
# are the only annotation that recognises the red-cell membrane skeleton. Sources are listed
# in PROVENANCE.md.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(purrr)
  library(readr)
  library(stringr)
})

uniprot_pattern <- "[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}"
first_acc <- \(x) sub("-[0-9]+$", "", str_extract(x, uniprot_pattern))
first_gene <- \(x) str_extract(x, "^[A-Za-z0-9][A-Za-z0-9-]*")

# CB2019 lists UniProt entry names (SLC4A1_HUMAN); the other four carry accession and gene.
parsers <- list(
  RESPIRE = \(d) transmute(d, acc = first_acc(Accession), gene = first_gene(Gene)),
  Uniprot = \(d) transmute(d, acc = first_acc(Entry), gene = first_gene(`Gene Names (primary)`)),
  JPR2017 = \(d) transmute(d, acc = first_acc(`Protein IDs`), gene = first_gene(`Gene names`)),
  CB2019 = \(d) transmute(d, acc = NA_character_, gene = sub("_HUMAN$", "", `Entry name`)),
  `This study` = \(d) transmute(d, acc = first_acc(Accession), gene = first_gene(Gene))
)

rbc_reference <- imap(parsers, \(parse, sheet) {
  readxl::read_excel(
    here("00_Input", "RBC_proteins.xlsx"),
    sheet = sheet, .name_repair = "unique_quiet"
  ) |>
    parse() |>
    mutate(source = sheet)
}) |>
  list_rbind() |>
  filter(!is.na(acc) | !is.na(gene))
# CB2019 has no accession. Borrow one from the other sheets by gene, or its rows never merge with
# the same protein listed elsewhere and n_sources undercounts.
gene_acc <- rbc_reference |>
  filter(!is.na(acc), !is.na(gene)) |>
  distinct(gene, acc_from_gene = acc) |>
  slice_head(n = 1, by = gene)
rbc_reference <- rbc_reference |>
  left_join(gene_acc, by = "gene") |>
  mutate(key = coalesce(acc, acc_from_gene, gene)) |>
  summarise(
    acc = first(na.omit(c(acc, acc_from_gene))),
    gene = first(na.omit(gene)),
    sources = paste(sort(unique(source)), collapse = ";"),
    n_sources = n_distinct(source),
    .by = key
  ) |>
  select(acc, gene, sources, n_sources) |>
  arrange(desc(n_sources), coalesce(acc, gene))

write_tsv(rbc_reference, here("00_Input", "RBC_proteome_reference.tsv"))
message(
  "RBC_proteome_reference.tsv: ", nrow(rbc_reference), " proteins, ",
  sum(rbc_reference$n_sources >= 2), " backed by two or more sources"
)

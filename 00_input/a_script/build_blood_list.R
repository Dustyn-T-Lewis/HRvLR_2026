# Collapses the five published red-cell proteome sources in RBC_proteins.xlsx into one lookup:
# one row per protein, the sources that list it, and how many agree. Mature erythrocytes are
# enucleate, so no transcriptomic atlas represents their proteome; these MS references are the
# only way to recognise the red-cell membrane skeleton by annotation. See 00_input/PROVENANCE.md.

read_rbc_sources <- function(xlsx_path) {
  uniprot_re <- "[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}"
  first_acc <- function(x) sub("-[0-9]+$", "", stringr::str_extract(x, uniprot_re))
  first_gene <- function(x) stringr::str_extract(x, "^[A-Za-z0-9][A-Za-z0-9-]*")

  # CB2019 lists UniProt entry names (SLC4A1_HUMAN); the others carry accession + gene columns.
  sheets <- list(
    RESPIRE      = \(d) dplyr::transmute(d, acc = first_acc(Accession), gene = first_gene(Gene)),
    Uniprot      = \(d) dplyr::transmute(d, acc = first_acc(Entry), gene = first_gene(`Gene Names (primary)`)),
    JPR2017      = \(d) dplyr::transmute(d, acc = first_acc(`Protein IDs`), gene = first_gene(`Gene names`)),
    CB2019       = \(d) dplyr::transmute(d, acc = NA_character_, gene = sub("_HUMAN$", "", `Entry name`)),
    `This study` = \(d) dplyr::transmute(d, acc = first_acc(Accession), gene = first_gene(Gene))
  )

  purrr::imap(sheets, \(parse, sheet) {
    d <- readxl::read_excel(xlsx_path, sheet = sheet, .name_repair = "unique_quiet")
    parse(d) |> dplyr::mutate(source = sheet)
  }) |>
    purrr::list_rbind() |>
    dplyr::filter(!is.na(acc) | !is.na(gene))
}

build_rbc_reference <- function(xlsx_path) {
  read_rbc_sources(xlsx_path) |>
    dplyr::mutate(key = dplyr::coalesce(acc, gene)) |>
    dplyr::summarise(
      acc = dplyr::first(stats::na.omit(acc)),
      gene = dplyr::first(stats::na.omit(gene)),
      sources = paste(sort(unique(source)), collapse = ";"),
      n_sources = dplyr::n_distinct(source),
      .by = key
    ) |>
    dplyr::select(acc, gene, sources, n_sources) |>
    dplyr::arrange(dplyr::desc(n_sources), dplyr::coalesce(acc, gene))
}

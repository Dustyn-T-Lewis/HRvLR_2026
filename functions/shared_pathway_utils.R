# Gene-set construction shared by the pathway and network stages.

pacman::p_load(here, dplyr, tibble)

# Minimum detected members a set needs. fgsea, fry, singscore and module ORA
# all apply it to detected members, so a 200-member set with 3 measured
# proteins never scores like a fully covered one.
SET_FLOOR <- 15L
SET_CEILING <- 500L

# The biological-process terms of the GO Consortium's generic slim, read from
# the tracked .obo rather than typed out, so the list moves with the release.
read_goslim_bp <- function(obo = here("00_input", "goslim_generic.obo")) {
  ids <- GSEABase::ids(GSEABase::getOBOCollection(obo))
  ids <- ids[AnnotationDbi::Ontology(ids) %in% "BP"]
  tibble(go_id = ids, name = unname(AnnotationDbi::Term(ids)))
}

# Each GO term's theme is its most specific slim ancestor: of the slim terms
# that contain it (itself included), the one with the fewest descendants.
# Terms under no slim term get NA.
goslim_theme <- function(go_ids, slim = read_goslim_bp()) {
  offspring <- AnnotationDbi::mget(
    slim$go_id, GO.db::GOBPOFFSPRING,
    ifnotfound = NA
  )
  members <- purrr::map2(slim$go_id, offspring, \(s, o) {
    tibble(slim_id = s, go_id = unique(c(s, stats::na.omit(o))))
  }) |>
    purrr::list_rbind() |>
    mutate(breadth = n(), .by = "slim_id")
  best <- members |>
    filter(.data$go_id %in% go_ids) |>
    slice_min(.data$breadth, n = 1, by = "go_id", with_ties = FALSE) |>
    left_join(slim, by = c(slim_id = "go_id"))
  best$name[match(go_ids, best$go_id)]
}

# One set per slim term: every gene annotated to it or any descendant.
goslim_sets <- function(slim = read_goslim_bp()) {
  hits <- AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = slim$go_id, keytype = "GOALL", columns = "SYMBOL"
  )
  sets <- split(hits$SYMBOL, hits$GOALL) |> lapply(unique)
  label <- slim$name[match(names(sets), slim$go_id)]
  names(sets) <- paste0("GOSLIM_", toupper(gsub("[^A-Za-z0-9]+", "_", label)))
  sets
}

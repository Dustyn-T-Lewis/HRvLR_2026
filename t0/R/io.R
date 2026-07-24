# Load and validate the three inputs: the normalized protein matrix, the per-sample design,
# and the per-subject phenotype that carries both prediction targets. Every path comes from config.

load_proteins <- function(cfg) {
  raw <- readr::read_csv(here::here(cfg$paths$proteins), show_col_types = FALSE)
  annot_cols <- c(cfg$ids$protein_key, "protein", cfg$ids$gene_key, "description")
  annot_cols <- intersect(annot_cols, names(raw))
  annotation <- raw[, annot_cols]
  expr <- as.matrix(raw[, setdiff(names(raw), annot_cols)])
  rownames(expr) <- raw[[cfg$ids$protein_key]]
  list(expr = expr, annotation = annotation)
}

load_meta <- function(cfg) {
  meta <- readr::read_csv(here::here(cfg$paths$metadata), show_col_types = FALSE) |>
    janitor::clean_names(case = "none")
  pheno <- readr::read_csv(here::here(cfg$paths$phenotype), show_col_types = FALSE)

  resp_src <- cfg$targets$responder$source
  if (!identical(cfg$targets$responder$split, "median")) {
    stop("responder split '", cfg$targets$responder$split, "' not supported; use 'median'")
  }
  cut <- stats::median(pheno[[resp_src]], na.rm = TRUE)
  pheno <- pheno |>
    dplyr::mutate(
      responder = factor(
        dplyr::if_else(.data[[resp_src]] >= cut,
          cfg$targets$responder$levels[[2]], cfg$targets$responder$levels[[1]]
        ),
        levels = cfg$targets$responder$levels
      ),
      supplement = factor(.data[[cfg$targets$supplement$source]])
    ) |>
    dplyr::select(subject, responder, supplement,
      comp_hypertrophy = dplyr::all_of(resp_src),
      dplyr::starts_with("d_")
    )

  meta |>
    dplyr::transmute(
      sample = .data[[cfg$ids$sample_key]],
      subject = .data[[cfg$ids$subject_key]],
      group = .data[[cfg$ids$group_col]],
      timepoint = .data[[cfg$ids$time_col]]
    ) |>
    dplyr::left_join(pheno, by = "subject")
}

# Hallmark gene sets as symbol lists, and the protein-id -> gene map the pathway/module steps use.
load_gene_sets <- function(cfg) {
  args <- list(species = "Homo sapiens")
  m <- tryCatch(
    do.call(msigdbr::msigdbr, c(args, list(collection = cfg$gene_sets$collection))),
    error = function(e) do.call(msigdbr::msigdbr, c(args, list(category = cfg$gene_sets$collection)))
  )
  sets <- split(m$gene_symbol, m$gs_name)
  sizes <- lengths(sets)
  sets[sizes >= cfg$gene_sets$min_size & sizes <= cfg$gene_sets$max_size]
}

gene_map <- function(proteins, cfg) {
  setNames(
    proteins$annotation[[cfg$ids$gene_key]],
    proteins$annotation[[cfg$ids$protein_key]]
  )
}

validate_inputs <- function(proteins, meta, cfg) {
  expr <- proteins$expr
  stopifnot(
    "every measured sample needs metadata" =
      all(colnames(expr) %in% meta$sample),
    "protein ids must be unique" = !anyDuplicated(rownames(expr)),
    "timepoints missing" =
      all(cfg$design$timepoints %in% meta$timepoint),
    "responder target has NA" = !anyNA(meta$responder),
    "supplement target has NA" = !anyNA(meta$supplement),
    "permutation reps must be a positive number" =
      is.numeric(cfg$permutation$reps) && cfg$permutation$reps >= 1
  )
  invisible(TRUE)
}

# Custom recipe steps that make pathway and module features fold-aware. Both consume the protein
# columns and replace them: step_singscore scores gene sets per subject (leakage-safe by design),
# step_wgcna_me learns modules on the training rows in prep() and projects new rows in bake(). Placing
# the network fit inside the recipe is what stops the whole-cohort WGCNA leakage seen in F05.

step_singscore <- function(recipe, ..., role = "predictor", trained = FALSE,
                           gene_sets = NULL, gene_map = NULL, min_size = 10,
                           inputs = NULL, skip = FALSE,
                           id = recipes::rand_id("singscore")) {
  recipes::add_step(recipe, step_singscore_new(
    terms = rlang::enquos(...), role = role, trained = trained, gene_sets = gene_sets,
    gene_map = gene_map, min_size = min_size, inputs = inputs, skip = skip, id = id
  ))
}

step_singscore_new <- function(terms, role, trained, gene_sets, gene_map, min_size,
                               inputs, skip, id) {
  recipes::step(
    subclass = "singscore", terms = terms, role = role, trained = trained,
    gene_sets = gene_sets, gene_map = gene_map, min_size = min_size,
    inputs = inputs, skip = skip, id = id
  )
}

prep.step_singscore <- function(x, training, info = NULL, ...) {
  cols <- recipes::recipes_eval_select(x$terms, training, info)
  present <- unique(x$gene_map[cols])
  usable <- x$gene_sets[
    vapply(x$gene_sets, \(g) length(intersect(g, present)) >= x$min_size, logical(1))
  ]
  step_singscore_new(
    terms = x$terms, role = x$role, trained = TRUE, gene_sets = usable,
    gene_map = x$gene_map, min_size = x$min_size, inputs = cols, skip = x$skip, id = x$id
  )
}

# Scores exactly the gene sets fixed in prep, so train and held-out rows share columns. Upstream
# median imputation removes NA, and the matrix cast keeps a single held-out subject two-dimensional.
bake.step_singscore <- function(object, new_data, ...) {
  cols <- object$inputs
  mat <- t(as.matrix(new_data[, cols, drop = FALSE]))
  rownames(mat) <- object$gene_map[cols]
  mat <- mat[!is.na(rownames(mat)) & !duplicated(rownames(mat)), , drop = FALSE]

  ranked <- singscore::rankGenes(mat)
  scores <- vapply(object$gene_sets, function(genes) {
    singscore::simpleScore(ranked, upSet = intersect(genes, rownames(mat)))$TotalScore
  }, numeric(ncol(mat)))
  scores <- matrix(scores,
    nrow = ncol(mat),
    dimnames = list(colnames(mat), names(object$gene_sets))
  )

  dplyr::bind_cols(
    new_data[, setdiff(names(new_data), cols), drop = FALSE],
    tibble::as_tibble(scores)
  )
}

step_wgcna_me <- function(recipe, ..., role = "predictor", trained = FALSE,
                          wgcna = NULL, loadings = NULL, inputs = NULL, skip = FALSE,
                          id = recipes::rand_id("wgcna_me")) {
  recipes::add_step(recipe, step_wgcna_me_new(
    terms = rlang::enquos(...), role = role, trained = trained, wgcna = wgcna,
    loadings = loadings, inputs = inputs, skip = skip, id = id
  ))
}

step_wgcna_me_new <- function(terms, role, trained, wgcna, loadings, inputs, skip, id) {
  recipes::step(
    subclass = "wgcna_me", terms = terms, role = role, trained = trained,
    wgcna = wgcna, loadings = loadings, inputs = inputs, skip = skip, id = id
  )
}

prep.step_wgcna_me <- function(x, training, info = NULL, ...) {
  cols <- recipes::recipes_eval_select(x$terms, training, info)
  frame <- tibble::as_tibble(training[, cols, drop = FALSE]) |>
    dplyr::mutate(subject = as.character(seq_len(nrow(training))), .before = 1)
  mods <- fit_wgcna_modules(frame, x$wgcna)
  step_wgcna_me_new(
    terms = x$terms, role = x$role, trained = TRUE, wgcna = x$wgcna,
    loadings = mods, inputs = cols, skip = x$skip, id = x$id
  )
}

bake.step_wgcna_me <- function(object, new_data, ...) {
  frame <- tibble::as_tibble(new_data[, object$inputs, drop = FALSE]) |>
    dplyr::mutate(subject = as.character(seq_len(nrow(new_data))), .before = 1)
  me <- project_eigengenes(object$loadings, frame) |> dplyr::select(-subject)
  dplyr::bind_cols(
    new_data[, setdiff(names(new_data), object$inputs), drop = FALSE], me
  )
}

print.step_singscore <- function(x, ...) {
  n <- if (x$trained) length(x$gene_sets) else NA_integer_
  cat("singscore gene-set scoring [", n, "sets ]\n")
  invisible(x)
}

print.step_wgcna_me <- function(x, ...) {
  n <- if (x$trained) length(x$loadings$loadings) else NA_integer_
  cat("WGCNA module eigengenes [", n, "modules ]\n")
  invisible(x)
}

tidy.step_singscore <- function(x, ...) {
  terms <- if (x$trained) names(x$gene_sets) else NA_character_
  tibble::tibble(terms = terms, id = x$id)
}

tidy.step_wgcna_me <- function(x, ...) {
  terms <- if (x$trained) paste0("ME", names(x$loadings$loadings)) else NA_character_
  tibble::tibble(terms = terms, id = x$id)
}

required_pkgs.step_singscore <- function(x, ...) "singscore"
required_pkgs.step_wgcna_me <- function(x, ...) c("WGCNA", "matrixStats")

# Register the S3 methods so recipes dispatches to them when sourced outside a package.
local({
  reg <- function(generic, class, method) {
    if (requireNamespace("recipes", quietly = TRUE)) {
      registerS3method(generic, class, method, envir = asNamespace("recipes"))
    }
  }
  reg("prep", "step_singscore", prep.step_singscore)
  reg("bake", "step_singscore", bake.step_singscore)
  reg("print", "step_singscore", print.step_singscore)
  reg("tidy", "step_singscore", tidy.step_singscore)
  reg("required_pkgs", "step_singscore", required_pkgs.step_singscore)
  reg("prep", "step_wgcna_me", prep.step_wgcna_me)
  reg("bake", "step_wgcna_me", bake.step_wgcna_me)
  reg("print", "step_wgcna_me", print.step_wgcna_me)
  reg("tidy", "step_wgcna_me", tidy.step_wgcna_me)
  reg("required_pkgs", "step_wgcna_me", required_pkgs.step_wgcna_me)
})

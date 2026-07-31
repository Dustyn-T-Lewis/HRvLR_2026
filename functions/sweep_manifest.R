# One row per swept cell: its metric at the largest B, at B = 0, its
# permutation p, whether it leads, and its named drivers. Derived from the store
# on demand rather than written to a file -- the store already holds every
# number here, so a saved manifest was a second copy of it keyed the same way.

pacman::p_load(here, dplyr)
source(here("functions", "sweep_grid.R"))

LEAKAGE_LABEL <- c(
  pathways = "leakage-free", modules = "optimistic", proteins = "optimistic"
)

METRIC_COLS <- list(
  class = c(metric = "estimate", p = "perm_p"),
  cont = c(metric = "q2", p = "perm_p_q2")
)

base_feature <- function(x) sub("@T[123]$", "", x)

manifest_drivers <- function(sel, n = 10) {
  if (is.null(sel) || !nrow(sel)) {
    return(list(n = 0L, top = ""))
  }
  pooled <- sel |>
    mutate(feature = base_feature(.data$feature)) |>
    group_by(.data$feature) |>
    summarise(freq = max(.data$freq), .groups = "drop") |>
    arrange(desc(.data$freq))
  list(
    n = nrow(pooled),
    top = paste(utils::head(pooled$feature, n), collapse = "; ")
  )
}

pred_row <- function(s, kind) {
  cols <- METRIC_COLS[[kind]]
  top_b <- s[s$B == max(s$B), , drop = FALSE][1, ]
  zero_b <- s[s$B == min(s$B), , drop = FALSE][1, ]
  data.frame(
    n = top_b$n, n_features = top_b$p,
    metric_name = cols[["metric"]],
    metric_b200 = top_b[[cols[["metric"]]]],
    metric_b0 = zero_b[[cols[["metric"]]]],
    perm_p = top_b[[cols[["p"]]]],
    null_mean = if ("null_q2_mean" %in% names(top_b)) {
      top_b$null_q2_mean
    } else {
      top_b$null_mean
    },
    null_sd = if ("null_q2_sd" %in% names(top_b)) {
      top_b$null_q2_sd
    } else {
      top_b$null_sd
    },
    stringsAsFactors = FALSE
  )
}

lead_flag <- function(kind, metric, p) {
  is_lead(metric, p, kind)
}

build_manifest <- function(root, kind, screen_size, root_name,
                           root_dir = sweep_root_dir(root)) {
  summary_store <- read_sweep_store(root, "summary", root_dir)
  selection_store <- read_sweep_store(root, "selection", root_dir)
  # No sort here: the store is written in key order already, and re-sorting with
  # dplyr would reorder it, because dplyr collates in the C locale and the store
  # follows the session locale the leaf-directory glob used.
  cells <- distinct(
    summary_store, .data$level, .data$config, .data$phenotype, .data$model
  )

  rows <- lapply(seq_len(nrow(cells)), function(i) {
    key <- cells[i, ]
    in_cell <- function(d) {
      if (is.null(d)) {
        return(NULL)
      }
      filter(
        d,
        .data$level == key$level, .data$config == key$config,
        .data$phenotype == key$phenotype, .data$model == key$model
      )
    }
    metrics <- pred_row(as.data.frame(in_cell(summary_store)), kind)
    dr <- manifest_drivers(in_cell(selection_store))

    cbind(
      as.data.frame(key, stringsAsFactors = FALSE),
      metrics,
      data.frame(
        is_lead = lead_flag(kind, metrics$metric_b200, metrics$perm_p),
        n_cells_screened = screen_size,
        leakage = LEAKAGE_LABEL[[key$level]],
        n_drivers = dr$n, top_drivers = dr$top,
        stringsAsFactors = FALSE
      )
    )
  })

  manifest <- bind_rows(rows)
  message(sprintf("%s manifest: %d cells", root_name, nrow(manifest)))
  manifest
}

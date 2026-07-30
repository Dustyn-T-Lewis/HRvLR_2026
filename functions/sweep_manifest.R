pacman::p_load(here, dplyr, openxlsx)
source(here("functions", "sweep_grid.R"))
source(here("functions", "sweep_cell_panel.R"))

LEAKAGE_LABEL <- c(
  pathways = "leakage-free", modules = "optimistic", proteins = "optimistic"
)

METRIC_COLS <- list(
  class = c(metric = "estimate", p = "perm_p"),
  cont = c(metric = "q2", p = "perm_p_q2")
)

base_feature <- function(x) sub("@T[123]$", "", x)

manifest_drivers <- function(path, n = 10) {
  if (!"selection" %in% openxlsx::getSheetNames(path)) {
    return(list(n = 0L, top = ""))
  }
  sel <- openxlsx::read.xlsx(path, "selection")
  if (!nrow(sel)) {
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

pred_row <- function(f, kind) {
  cols <- METRIC_COLS[[kind]]
  s <- openxlsx::read.xlsx(f, "summary")
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

build_manifest <- function(root, kind, screen_size, root_name) {
  files <- Sys.glob(
    file.path(root, "*", "*", "*", "*", "c_data", "results.xlsx")
  )

  rows <- lapply(files, function(f) {
    parts <- strsplit(f, .Platform$file.sep)[[1]]
    n <- length(parts)
    level <- parts[[n - 5]]
    metrics <- pred_row(f, kind)
    dr <- manifest_drivers(f)

    cbind(
      data.frame(
        level = level, config = parts[[n - 4]],
        phenotype = parts[[n - 3]], model = parts[[n - 2]],
        stringsAsFactors = FALSE
      ),
      metrics,
      data.frame(
        is_lead = lead_flag(kind, metrics$metric_b200, metrics$perm_p),
        n_cells_screened = screen_size,
        leakage = LEAKAGE_LABEL[[level]],
        n_drivers = dr$n, top_drivers = dr$top,
        panel_path = file.path(dirname(dirname(f)), "b_reports", "panel.png"),
        data_path = f,
        stringsAsFactors = FALSE
      )
    )
  })

  manifest <- bind_rows(rows)
  write_sweep_workbook(
    file.path(root, "MANIFEST.xlsx"), list(manifest = manifest)
  )
  message(sprintf("%s manifest: %d cells", root_name, nrow(manifest)))
  manifest
}

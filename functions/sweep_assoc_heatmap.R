# Shared row/tile plumbing for the association heatmaps. The F04 per-method
# grids and the F07 synthesis figure read the same leaf workbooks through these,
# so a change to row selection or the BH cut cannot land in one and not the
# other.

pacman::p_load(here, dplyr, openxlsx)
source(here("functions", "sweep_grid.R"))

BH_CUT <- 0.05
N_SHOW <- c(modules = Inf, pathways = 18, proteins = 18)

OUTCOME_SHORT <- c(
  comp_hypertrophy = "comp.hyp", d_fcsa_I = "fCSA I", d_fcsa_II = "fCSA II",
  d_mcsa = "mCSA", d_1rm_legpress = "1RM leg", d_1rm_ext = "1RM ext",
  group_diff = "HR-LR"
)

# One row per feature x outcome x config. The trajectory config stores each
# feature three times as feature@T1/@T2/@T3. Collapsing keeps the best of the
# three, which is a maximum, so `collapse_tp = FALSE` returns all three copies
# and lets a caller give each timepoint its own column instead.
assoc_grid <- function(level, method, collapse_tp = TRUE) {
  files <- Sys.glob(file.path(
    sweep_root_dir("F04_association"), level, "*", "*", method,
    "c_data", "results.xlsx"
  ))
  bind_rows(lapply(files, function(f) {
    read.xlsx(f, "cell") |>
      transmute(
        tp = ifelse(
          grepl("@T[123]$", .data$feature), sub("^.*@", "", .data$feature), NA
        ),
        feature = sub("@T[123]$", "", .data$feature),
        n = .data$n, effect = .data$effect, t = .data$t, p = .data$p,
        bh = .data$bh,
        outcome = basename(dirname(dirname(dirname(f)))),
        config = basename(dirname(dirname(dirname(dirname(f)))))
      )
  })) |>
    (\(d) if (collapse_tp) collapse_timepoints(d) else d)()
}

collapse_timepoints <- function(d) {
  d |>
    group_by(.data$feature, .data$outcome, .data$config) |>
    slice_min(.data$p, n = 1, with_ties = FALSE) |>
    ungroup()
}

# Every BH-surviving feature earns a row, then the strongest nominal features
# fill the rest. Ranking on p alone can drop a survivor, because BH is applied
# within a cell and the cells hold 11 to 5,711 features.
heatmap_features <- function(d, level) {
  n_show <- N_SHOW[[level]]
  if (is.infinite(n_show)) {
    return(unique(d$feature))
  }
  survivors <- unique(d$feature[d$bh < BH_CUT])
  ranked <- d |>
    filter(!.data$feature %in% survivors) |>
    group_by(.data$feature) |>
    summarise(best = min(.data$p, na.rm = TRUE), .groups = "drop") |>
    slice_min(.data$best,
      n = max(0, n_show - length(survivors)),
      with_ties = FALSE
    ) |>
    pull(.data$feature)
  c(survivors, ranked)
}

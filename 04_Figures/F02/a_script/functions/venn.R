# Baseline-anchored 3-set Euler/venn helpers shared by F02 panels I and J.
# Globals resolved at call time from HRvLR_F02_setup.R: dep_df, PANEL_H, RPT_DIR,
# FIG_THEME, FIG_GEOM_TEXT, FIG_AXIS_TEXT, DIR_COLORS, GROUP_COLORS, CONTRAST_COLORS, save_png.

pacman::p_load(here, dplyr, tidyr, tibble, ggplot2, ggforce, patchwork, eulerr)

in_ellipse <- function(x, y, e) {
  dx <- x - e$h
  dy <- y - e$k
  xr <- dx * cos(e$phi) + dy * sin(e$phi)
  yr <- -dx * sin(e$phi) + dy * cos(e$phi)
  (xr / e$a)^2 + (yr / e$b)^2 <= 1
}

# Centroid of each disjoint region, found by classifying a fine grid over the
# ellipse bounding box; falls back to the mean of member-set centres if a region
# has no grid coverage (high-stress fits).
region_centroids <- function(ell, set_names) {
  pad <- 0.15 * max(ell$a, ell$b)
  xs <- seq(min(ell$h - ell$a) - pad, max(ell$h + ell$a) + pad, length.out = 220)
  ys <- seq(min(ell$k - ell$b) - pad, max(ell$k + ell$b) + pad, length.out = 220)
  grid <- expand.grid(x = xs, y = ys)
  member <- vapply(set_names, function(s) {
    e <- as.list(ell[s, ])
    in_ellipse(grid$x, grid$y, e)
  }, logical(nrow(grid)))
  key <- apply(member, 1, function(m) paste(set_names[m], collapse = "&"))
  grid$key <- key
  cent <- grid |>
    filter(key != "") |>
    group_by(key) |>
    summarise(cx = mean(x), cy = mean(y), .groups = "drop")
  set_center <- function(k) {
    s <- strsplit(k, "&")[[1]]
    c(cx = mean(ell[s, "h"]), cy = mean(ell[s, "k"]))
  }
  list(grid = cent, fallback = set_center)
}

# 3-set (or 2-set) Euler drawn in ggplot with per-region counts; no gene labels.
euler_only <- function(sets, cols, title, subtitle) {
  set_names <- names(sets)
  fit <- euler(sets, shape = "ellipse")
  ell <- as.data.frame(fit$ellipses)
  ell$set <- rownames(ell)

  all_genes <- unique(unlist(sets))
  membership <- vapply(
    set_names, function(s) all_genes %in% sets[[s]],
    logical(length(all_genes))
  )
  region_key <- apply(membership, 1, function(m) paste(set_names[m], collapse = "&"))
  gene_df <- tibble(gene = all_genes, key = region_key) |> filter(key != "")

  cents <- region_centroids(ell, set_names)
  centroid_xy <- function(k) {
    hit <- cents$grid |> filter(key == k)
    if (nrow(hit)) c(cx = hit$cx, cy = hit$cy) else cents$fallback(k)
  }
  region_counts <- gene_df |>
    count(key, name = "n") |>
    rowwise() |>
    mutate(cx = centroid_xy(key)[["cx"]], cy = centroid_xy(key)[["cy"]]) |>
    ungroup()

  gx <- mean(ell$h)
  gy <- mean(ell$k)
  name_df <- ell |>
    mutate(
      dx = h - gx, dy = k - gy,
      norm = pmax(sqrt(dx^2 + dy^2), 1e-6),
      reach = pmax(a, b) * 1.15,
      lx = h + dx / norm * reach,
      ly = k + dy / norm * reach
    )

  ggplot() +
    ggforce::geom_ellipse(
      data = ell,
      aes(x0 = h, y0 = k, a = a, b = b, angle = phi, fill = set, color = set),
      alpha = 0.32, linewidth = 0.7
    ) +
    geom_label(
      data = region_counts, aes(cx, cy, label = n),
      size = FIG_GEOM_TEXT, fontface = "bold", color = "grey15",
      fill = scales::alpha("white", 0.75), label.size = 0,
      label.padding = unit(0.8, "pt")
    ) +
    geom_text(
      data = name_df, aes(lx, ly, label = set, color = set),
      fontface = "bold", size = FIG_GEOM_TEXT, show.legend = FALSE
    ) +
    scale_fill_manual(values = cols, guide = "none") +
    scale_color_manual(values = cols, guide = "none") +
    coord_equal(clip = "off") +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    FIG_THEME +
    theme(
      axis.text = element_blank(), axis.ticks = element_blank(),
      panel.grid = element_blank(), panel.border = element_blank(),
      plot.margin = margin(4, 2, 4, 6)
    )
}

# Companion diverging bar: down (left, blue) / up (right, red) protein counts per
# overlap region. Direction is the sign of logFC in the region's defining contrast.
direction_bar <- function(bar_df) {
  lim <- max(bar_df$n, 1) * 1.4
  ggplot(bar_df, aes(signed, region, fill = direction)) +
    geom_col(width = 0.66, color = "white", linewidth = 0.25) +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "grey40") +
    geom_text(
      data = filter(bar_df, n > 0),
      aes(label = n, hjust = if_else(direction == "Down", 1.25, -0.25)),
      size = FIG_GEOM_TEXT, fontface = "bold", color = "grey20"
    ) +
    scale_fill_manual(
      values = c(Down = unname(DIR_COLORS[["Down"]]), Up = unname(DIR_COLORS[["Up"]])),
      breaks = c("Down", "Up"), name = NULL
    ) +
    scale_x_continuous(
      limits = c(-lim, lim), labels = abs,
      breaks = scales::pretty_breaks(3)
    ) +
    labs(
      title = "Direction", subtitle = "down / up by region contrast",
      x = "proteins", y = NULL
    ) +
    FIG_THEME +
    theme(
      legend.position = "bottom",
      axis.text.y = element_text(face = "bold", size = FIG_AXIS_TEXT - 1),
      panel.grid.major.y = element_blank(),
      plot.margin = margin(4, 6, 4, 2)
    )
}

# Euler + companion bar stitched into one saved panel.
euler_gg <- function(sets, cols, bar_df, title, subtitle, file_stem, w = 150, h = PANEL_H) {
  combined <- euler_only(sets, cols, title, subtitle) +
    direction_bar(bar_df) +
    plot_layout(widths = c(0.64, 0.36))
  save_png(combined, file_stem, w, h)
  invisible(combined)
}

# Genes passing a metric for one contrast.
dep_genes <- function(contrast, metric = c("pi", "p")) {
  metric <- match.arg(metric)
  sel <- if (metric == "pi") {
    dep_df[[paste0("sig_pi_", contrast)]] != 0
  } else {
    dep_df[[paste0("P.Value_", contrast)]] < 0.05
  }
  keep <- which(sel & !is.na(dep_df$gene))
  unique(dep_df$gene[keep])
}

# Disjoint up/down counts for the four companion-bar regions. Direction uses the
# region's defining contrast: shared/HR-only -> HR, LR-only -> LR, Baseline -> base.
region_bar_df <- function(hr_ct, lr_ct, base_ct, metric) {
  sel <- function(ct) {
    if (metric == "pi") {
      dep_df[[paste0("sig_pi_", ct)]] != 0
    } else {
      dep_df[[paste0("P.Value_", ct)]] < 0.05
    }
  }
  regions <- c("Shared", "HR-only", "LR-only", "Baseline")
  tibble(
    gene     = dep_df$gene,
    in_hr    = sel(hr_ct) %in% TRUE,
    in_lr    = sel(lr_ct) %in% TRUE,
    in_base  = sel(base_ct) %in% TRUE,
    lfc_hr   = dep_df[[paste0("logFC_", hr_ct)]],
    lfc_lr   = dep_df[[paste0("logFC_", lr_ct)]],
    lfc_base = dep_df[[paste0("logFC_", base_ct)]]
  ) |>
    filter(!is.na(gene), in_hr | in_lr | in_base) |>
    distinct(gene, .keep_all = TRUE) |>
    mutate(
      region = case_when(
        in_hr & in_lr ~ "Shared",
        in_hr & !in_lr ~ "HR-only",
        in_lr & !in_hr ~ "LR-only",
        in_base ~ "Baseline"
      ),
      lfc_use = case_when(
        region %in% c("Shared", "HR-only") ~ lfc_hr,
        region == "LR-only" ~ lfc_lr,
        region == "Baseline" ~ lfc_base
      ),
      direction = if_else(lfc_use >= 0, "Up", "Down")
    ) |>
    filter(!is.na(region), !is.na(lfc_use)) |>
    count(region, direction, name = "n") |>
    complete(region = regions, direction = c("Down", "Up"), fill = list(n = 0)) |>
    mutate(
      signed = if_else(direction == "Down", -n, n),
      region = factor(region, levels = rev(regions))
    )
}

venn3_panel <- function(hr_ct, lr_ct, base_ct, title, file_stem,
                        metric = c("pi", "p"), w = 150, h = PANEL_H) {
  metric <- match.arg(metric)
  sets <- list(
    HR = dep_genes(hr_ct, metric),
    LR = dep_genes(lr_ct, metric),
    Baseline = dep_genes(base_ct, metric)
  )
  cols <- c(
    HR = unname(GROUP_COLORS[["HR"]]),
    LR = unname(GROUP_COLORS[["LR"]]),
    Baseline = unname(CONTRAST_COLORS[["Baseline_HRvLR"]])
  )
  metric_lab <- if (metric == "pi") "Pi < 0.05" else "p < 0.05"
  euler_gg(sets, cols, region_bar_df(hr_ct, lr_ct, base_ct, metric),
    title, metric_lab, file_stem,
    w = w, h = h
  )
}

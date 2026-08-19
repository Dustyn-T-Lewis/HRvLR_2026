# F02 supplement: pi-selected protein heatmaps, real arm labels beside shuffled ones
# One row per timepoint. Left column selects proteins by pi < PI_THRESH on the real arm
# contrast; right column runs the identical selection on labels shuffled across subjects.
# Both columns block cleanly, because selecting proteins for separating the arms and then
# displaying that separation is one fact shown twice. The shuffled column is the control
# that says so. Rows are grouped by GO Slim, the compact axis, so blocks can be read
# against biology rather than against the selection.

pacman::p_load(here, dplyr, tidyr, tibble, ggplot2, patchwork, limma, withr)

if (!exists("meta")) source(here("04_Figures", "F02_proteome", "a_script", "setup.R"))
source(here("03_Features", "contrasts.R"))
source(here("functions", "shared_pathway_utils.R"))

PI_HM_SEED <- 7
TP_LABEL <- c(T1 = "T1 baseline", T2 = "T2 trained", T3 = "T3 acute")

goslim <- build_goslim_gene_sets(min_size = SET_FLOOR, max_size = 500)
# A protein sits in several slim terms; take the smallest containing term so the label
# is the most specific one available, and leave the rest unassigned.
slim_of <- function(genes) {
  sizes <- lengths(goslim)
  vapply(genes, function(g) {
    hit <- names(goslim)[vapply(goslim, function(s) g %in% s, logical(1))]
    if (!length(hit)) "Unassigned" else hit[which.min(sizes[hit])]
  }, character(1))
}

pi_selected <- function(x, g) {
  tt <- topTable(eBayes(lmFit(x, model.matrix(~g))),
    coef = 2, number = Inf, sort.by = "none"
  )
  rownames(x)[pi_score(tt$P.Value, tt$logFC) < PI_THRESH]
}

heat_frame <- function(tp, labels, arm, tag) {
  samples <- meta$Col_ID[meta$Timepoint == tp]
  x <- imp_mat[, samples]
  keep <- pi_selected(x, factor(labels))
  z <- t(scale(t(x[keep, , drop = FALSE])))
  as.data.frame(z) |>
    rownames_to_column("gene") |>
    pivot_longer(-"gene", names_to = "sample", values_to = "z") |>
    mutate(
      slim = slim_of(.data$gene)[.data$gene],
      arm = arm[.data$sample],
      timepoint = TP_LABEL[[tp]],
      tag = tag,
      n_sel = length(keep)
    )
}

set_labels <- function(tp) {
  samples <- meta$Col_ID[meta$Timepoint == tp]
  real <- as.character(meta$Group[match(samples, meta$Col_ID)])
  list(
    samples = samples, real = real,
    shuffled = with_seed(PI_HM_SEED + match(tp, names(TP_LABEL)), sample(real))
  )
}

frames <- lapply(names(TP_LABEL), function(tp) {
  lab <- set_labels(tp)
  arm_real <- stats::setNames(lab$real, lab$samples)
  arm_shuf <- stats::setNames(lab$shuffled, lab$samples)
  bind_rows(
    heat_frame(tp, lab$real, arm_real, "Selected on arm"),
    heat_frame(tp, lab$shuffled, arm_shuf, "Selected on shuffle")
  )
})
pi_heat <- bind_rows(frames) |>
  mutate(
    timepoint = factor(.data$timepoint, levels = unname(TP_LABEL)),
    tag = factor(.data$tag, levels = c("Selected on arm", "Selected on shuffle"))
  )

# Order samples by the label the selection used, and proteins by GO Slim block.
pi_heat <- pi_heat |>
  arrange(.data$timepoint, .data$tag, .data$arm, .data$sample) |>
  group_by(.data$timepoint, .data$tag) |>
  mutate(
    sample = factor(.data$sample, levels = unique(.data$sample)),
    gene = factor(.data$gene, levels = unique(.data$gene[order(.data$slim, .data$gene)]))
  ) |>
  ungroup()

counts <- distinct(pi_heat, timepoint, tag, n_sel)

# Each cell gets its own axes. facet_grid would spread every panel's tiles over the
# union of all gene and sample levels, leaving six sparse rectangles.
one_heat <- function(tp, tg) {
  d <- filter(pi_heat, .data$timepoint == tp, .data$tag == tg) |>
    mutate(
      gene = factor(.data$gene, levels = unique(.data$gene[order(.data$slim, .data$gene)])),
      sample = factor(.data$sample, levels = unique(.data$sample[order(.data$arm, .data$sample)]))
    )
  slim_runs <- d |>
    distinct(.data$gene, .data$slim) |>
    arrange(.data$gene) |>
    count(.data$slim, name = "n")

  ggplot(d, aes(sample, gene, fill = z)) +
    geom_raster() +
    scale_fill_gradient2(
      low = DIR_COLORS[["Down"]], mid = "white", high = DIR_COLORS[["Up"]],
      midpoint = 0, limits = c(-2, 2), oob = scales::squish, name = "z"
    ) +
    labs(
      title = sprintf(
        "%s  |  %s  |  %d proteins, %d GO Slim groups",
        tp, tg, d$n_sel[1], nrow(slim_runs)
      ),
      x = NULL, y = NULL
    ) +
    FIG_THEME +
    theme(
      axis.text = element_blank(), axis.ticks = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(size = FIG_GEOM_TEXT, face = "bold", colour = "grey20")
    )
}

grid <- lapply(levels(pi_heat$timepoint), function(tp) {
  lapply(levels(pi_heat$tag), function(tg) one_heat(tp, tg))
})

p_pi <- wrap_plots(unlist(grid, recursive = FALSE), ncol = 2, guides = "collect") +
  plot_annotation(
    caption = paste(
      "Left selects proteins by pi < 0.05 on the real arm contrast; right runs the",
      "identical selection on labels shuffled across subjects.\nBoth block, and at every",
      "timepoint the shuffle selects more proteins than the arm does. Samples ordered by",
      "the label the selection used;\nrows grouped by GO Slim."
    ),
    theme = theme(
      plot.caption = element_text(hjust = 0, size = FIG_GEOM_TEXT - 0.4, colour = "grey35")
    )
  )

save_png(p_pi, file.path(RPT_DIR, "supp", "supp_pi_heatmap"), 200, 210)
F02_AUDIT[["supp_pi_heatmap"]] <- pi_heat |>
  distinct(timepoint, tag, gene, slim, n_sel)
cat("F02 pi heatmap supplement done.\n")

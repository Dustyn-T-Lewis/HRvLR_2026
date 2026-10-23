# The protein dot-matrix (F2): one row per protein, labelled by gene, and one column per contrast
# under its window strip. Colour is logFC (PiYG), size -log10 p; a filled dot is a Π hit
# (p^|logFC| < alpha_nominal), an open one is not, and a black outline marks FDR < alpha_fdr.
# A filled dot carries logFC as its fill, an open one as its outline. Colour is capped (see below);
# the logFC in the panel data is not.
# Each row reads "gene (observed/n_samples)", the biopsies where a precursor was seen. A protein
# seen in few biopsies takes most of its values from limpa's detection model, which can give it a
# large logFC on little data, and so a small Π; rows are therefore ranked by p, which carries
# limpa's quantification uncertainty, and Π only decides which proteins may enter.
# Expects dplyr and ggplot2 attached and figure_style.R sourced.
pi_calls <- c(sprintf("Π < %s", alpha_nominal), sprintf("Π ≥ %s", alpha_nominal))
fdr_call <- sprintf("FDR < %s", alpha_fdr)

# A few proteins (blood neutrophil markers after the bout) move far more than the rest; the colour
# stops at the median |logFC| of all Π hits, rounded up to a quarter, so the rest stay readable.
# Taken over the whole table so every protein panel shares one scale.
protein_colour_limit <- function(results) {
  ceiling(median(abs(results$logFC[results$pi_score < alpha_nominal])) * 4) / 4
}

# The HR - LR candidates: proteins with a Π hit in a between-group or interaction contrast, each
# placed in the block of its hit with the smallest p and ranked on that p.
hr_lr_candidates <- function(results) {
  results |>
    filter(
      contrast_panel(contrast) %in% c("Between groups", "Interaction"),
      pi_score < alpha_nominal
    ) |>
    slice_min(p, n = 1, by = uniprot_id, with_ties = FALSE) |>
    transmute(
      uniprot_id, gene, n_observed,
      block = droplevels(contrast_panel(contrast)), selected_in = contrast, selection_p = p,
      logFC, p, fdr, pi_score
    )
}

# Every contrast of the chosen proteins, the per_row with the smallest selection_p in each block,
# checked against the source table. `chosen` holds uniprot_id, block, selection_p and any extra
# columns to carry (such as selected_in).
protein_panel_data <- function(results, chosen, per_row) {
  contrasts <- unique(results$contrast)
  shown <- chosen |>
    slice_min(selection_p, n = per_row, by = block, with_ties = FALSE) |>
    mutate(row_order = rank(selection_p, ties.method = "first"), .by = block) |>
    select(uniprot_id, block, row_order, any_of("selected_in"), selection_p)
  panel_data <- results |>
    inner_join(shown, by = "uniprot_id") |>
    select(
      block, row_order, uniprot_id, gene, n_observed, contrast, logFC, logFC_lower, logFC_upper,
      p, fdr, pi_score, any_of("selected_in"), selection_p
    ) |>
    arrange(block, row_order, match(contrast, contrasts))
  check <- inner_join(panel_data, results, by = c("uniprot_id", "contrast"))
  stopifnot(
    nrow(panel_data) == length(contrasts) * nrow(shown), nrow(check) == nrow(panel_data),
    all(check$logFC.x == check$logFC.y), all(check$p.x == check$p.y),
    all(check$fdr.x == check$fdr.y), all(check$pi_score.x == check$pi_score.y),
    max(count(shown, block)$n) <= per_row
  )
  panel_data
}

# `results` has one row per chosen protein and contrast, with the protein_by_contrast columns plus
# `block` (row facet, a factor) and `row_order` (row position within the block, 1 at the top).
# `contrasts` gives the column order; `limit` the colour cap.
protein_grid <- function(results, contrasts, title, subtitle, n_samples, limit, p_cap = 4) {
  d <- results |>
    mutate(
      family = contrast_panel(contrast),
      column = factor(contrast, contrasts),
      row = reorder(paste(block, uniprot_id), -row_order),
      log_p = pmin(-log10(p), p_cap),
      pi_call = factor(if_else(pi_score < alpha_nominal, pi_calls[1], pi_calls[2]), pi_calls)
    )
  # A gene shared by two accessions gets the accession, so no two rows read alike.
  labels <- distinct(d, row, gene, uniprot_id, n_observed) |>
    mutate(
      accession = if_else(
        duplicated(gene) | duplicated(gene, fromLast = TRUE), paste0(uniprot_id, ", "), ""
      ),
      label = sprintf("%s (%s%d/%d)", gene, accession, as.integer(n_observed), n_samples)
    )
  strips <- contrast_strips(x = levels(d$family), y = levels(d$block))
  ggplot(d, aes(column, row)) +
    geom_point(
      aes(size = log_p, colour = logFC, shape = pi_call), filter(d, pi_call == pi_calls[2]),
      stroke = 0.4
    ) +
    geom_point(
      aes(size = log_p, fill = logFC, shape = pi_call), filter(d, pi_call == pi_calls[1]),
      colour = "grey30", stroke = 0.25
    ) +
    geom_point(
      aes(size = log_p, alpha = fdr_call), filter(d, fdr < alpha_fdr),
      shape = 1, colour = "black", stroke = 0.7
    ) +
    scale_colour_distiller(
      palette = "PiYG", direction = 1, limits = c(-limit, limit), oob = scales::squish,
      aesthetics = c("colour", "fill"),
      name = sprintf("logFC\n(capped at\n±%s)", limit)
    ) +
    scale_size_area(
      max_size = 2.6, limits = c(0, p_cap), breaks = c(1, 2, p_cap),
      labels = c("1", "2", paste0(p_cap, "+")), name = "-log10 p"
    ) +
    scale_shape_manual(values = rlang::set_names(c(21, 1), pi_calls), name = NULL, drop = FALSE) +
    scale_alpha_manual(values = rlang::set_names(1, fdr_call), name = NULL) +
    scale_x_discrete(labels = rlang::set_names(contrast_markup(contrasts), contrasts)) +
    scale_y_discrete(labels = tibble::deframe(select(labels, row, label))) +
    guides(
      colour = guide_colourbar(order = 1, barheight = unit(14, "mm"), barwidth = unit(2.5, "mm")),
      fill = "none",
      size = guide_legend(
        order = 2,
        override.aes = list(shape = 21, fill = "grey60", colour = "grey30", stroke = 0.25)
      ), shape = guide_legend(
        order = 3, override.aes = list(size = 2, fill = "grey60", colour = "grey30")
      ),
      alpha = guide_legend(order = 4, override.aes = list(size = 2))
    ) +
    ggh4x::facet_grid2(
      block ~ family,
      scales = "free", space = "free", switch = "y", strip = strips,
      labeller = labeller(family = label_wrap_gen(10))
    ) +
    labs(x = NULL, y = NULL, title = title, subtitle = subtitle) +
    figure_theme +
    matrix_theme +
    theme(axis.text.y = element_text(size = 5.5))
}

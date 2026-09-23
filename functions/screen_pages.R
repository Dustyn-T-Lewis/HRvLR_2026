# The two screen pages every feature level's packet carries, drawn the same
# way so proteins, pathways and modules read side by side.

pacman::p_load(here, dplyr, purrr, ggplot2)

source(here("functions", "classify.R"))
source(here("functions", "shared_style.R"))

auc_page <- function(classify, level, data_note) {
  classify |>
    filter(!is.na(.data$auc)) |>
    mutate(task = factor(.data$task, levels = TASKS$task)) |>
    ggplot(aes(.data$auc)) +
    geom_histogram(
      breaks = seq(0, 1, 0.05), fill = "grey55", colour = "white"
    ) +
    geom_vline(xintercept = 0.5, linetype = "dashed") +
    facet_wrap(~task, nrow = 2) +
    labs(
      title = sprintf("%s AUC per classification task", level),
      subtitle = sprintf(
        "Mann-Whitney AUC, Wilcoxon p (signed-rank within arm); %d features",
        n_distinct(classify$feature)
      ),
      x = "AUC (above 0.5 = higher at the later timepoint, or in HR)",
      y = "Features",
      caption = caption(
        "Each feature contributes one AUC per task. Within-arm tasks compare ",
        "a subject's later timepoint with its earlier one; HRvLR tasks ",
        "compare subjects on the baseline level or on their own change. ",
        "Dashed line: ",
        "chance. Data: ", data_note, ", sheets classify and chance_classify."
      )
    ) +
    FIG_THEME
}

association_page <- function(chance, level, data_note) {
  ggplot(chance, aes(.data$phenotype, .data$window, fill = .data$ratio)) +
    geom_tile(colour = "white") +
    geom_text(
      aes(label = sprintf("%d\n%d", .data$n_nominal, .data$n_bh)),
      size = 2.6
    ) +
    scale_fill_gradient2(
      low = "#4393C3", mid = "white", high = "#D6604D", midpoint = 1
    ) +
    scale_y_discrete(limits = c("acute", "training", "T1")) +
    labs(
      title = sprintf("%s association with phenotype", level),
      subtitle = sprintf(
        "limma on a continuous predictor; BH per cell; %d features",
        max(chance$n_tested)
      ),
      x = NULL, y = "Window", fill = "Nominal /\nexpected",
      caption = caption(
        "Rows: the feature at T1, or its change over training (T2 - T1) or ",
        "the acute bout (T3 - T2). Columns: the ten phenotypes. Fill: nominal ",
        "hits at p < 0.05 over the count expected by chance (0.05 x tested); ",
        "white is what an empty screen looks like. Top number: nominal hits; ",
        "bottom: BH < 0.05. Data: ", data_note,
        ", sheets associate and chance_associate."
      )
    ) +
    FIG_THEME +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
}

# Every feature that reaches nominal p < 0.05 in at least one column, as a dot
# matrix: rows are features, columns are tasks, contrasts or phenotypes. Rows
# are ordered by how many columns they are nominal in, then by their best p,
# and split into pages of `per_page`, so every nominal feature is drawn.
hit_pages <- function(res, col_levels, title, subtitle, effect_label,
                      data_note, midpoint = 0, per_page = 75L) {
  ranked <- res |>
    filter(!is.na(.data$p)) |>
    summarise(
      n_nominal = sum(.data$p < 0.05), best = min(.data$p),
      .by = "label"
    ) |>
    filter(.data$n_nominal > 0) |>
    arrange(desc(.data$n_nominal), .data$best)
  if (!nrow(ranked)) {
    return(list())
  }
  page_of <- ceiling(seq_len(nrow(ranked)) / per_page)
  n_pages <- max(page_of)
  pages <- map(seq_len(n_pages), function(pg) {
    rows <- ranked$label[page_of == pg]
    res |>
      filter(.data$label %in% rows) |>
      mutate(
        label = factor(.data$label, levels = rev(rows)),
        column = factor(.data$column, levels = col_levels),
        nominal = !is.na(.data$p) & .data$p < 0.05
      ) |>
      ggplot(aes(.data$column, .data$label)) +
      geom_point(
        data = \(d) filter(d, !.data$nominal),
        colour = "grey85", size = 0.6
      ) +
      geom_point(
        data = \(d) filter(d, .data$nominal),
        aes(size = -log10(.data$p), fill = .data$effect),
        shape = 21, colour = "grey30", stroke = 0.2
      ) +
      geom_point(
        data = \(d) filter(d, .data$bh < 0.05),
        aes(size = -log10(.data$p)),
        shape = 21, colour = "black", stroke = 0.9, fill = NA
      ) +
      scale_fill_gradient2(
        low = "#2166AC", mid = "white", high = "#B2182B",
        midpoint = midpoint
      ) +
      scale_size_continuous(range = c(1, 4.5), limits = c(-log10(0.05), NA)) +
      scale_x_discrete(drop = FALSE) +
      labs(
        title = sprintf("%s (%d/%d)", title, pg, n_pages),
        subtitle = sprintf(
          "%s; %d features nominal in at least one column; rows %d-%d",
          subtitle, nrow(ranked), min(which(page_of == pg)),
          max(which(page_of == pg))
        ),
        x = NULL, y = NULL, fill = effect_label, size = "-log10 p",
        caption = caption(
          "Rows: features at nominal p < 0.05 in at least one column, most ",
          "columns first, then by best p. Filled dot: nominal in that column, ",
          "fill = ", effect_label, ", size = -log10 p. Black ring: BH < 0.05. ",
          "Grey speck: tested, not nominal. ",
          "Data: ", data_note, "."
        )
      ) +
      FIG_THEME +
      theme(
        axis.text.x = element_text(angle = 35, hjust = 1),
        axis.text.y = element_text(size = if (length(rows) > 30) 5.5 else 8),
        panel.grid.major = element_line(colour = "grey94", linewidth = 0.2)
      )
  })
  set_names(pages, sprintf("%d/%d", seq_len(n_pages), n_pages))
}

# The hit pages a screen produces at any level, classification then
# association for each window, flattened to one named list of pages.
screen_hit_pages <- function(screens, level, data_note, label_of = identity) {
  classify <- screens$classify |>
    transmute(
      label = label_of(.data$feature), column = .data$task,
      effect = .data$auc, p = .data$p, bh = .data$bh
    )
  pages <- list(hit_pages(
    classify, TASKS$task,
    title = sprintf("%s classification hits", level),
    subtitle = "AUC and Wilcoxon p per task, paired within arm",
    effect_label = "AUC", midpoint = 0.5,
    data_note = paste0(data_note, ", sheet classify")
  ))
  names(pages) <- sprintf("%s classification hits", level)
  window_names <- c(
    T1 = "level at T1", training = "training change (T2 - T1)",
    acute = "acute change (T3 - T2)"
  )
  for (w in names(window_names)) {
    assoc <- screens$associate |>
      filter(.data$window == w) |>
      transmute(
        label = label_of(.data$feature), column = .data$phenotype,
        effect = .data$t, p = .data$p, bh = .data$bh
      )
    key <- sprintf("%s association hits, %s", level, window_names[[w]])
    pages[[key]] <- hit_pages(
      assoc, PHENOTYPES,
      title = key,
      subtitle = "limma slope t against each phenotype, one row per subject",
      effect_label = "t",
      data_note = paste0(data_note, ", sheet associate")
    )
  }
  list_flatten(pages, name_spec = "{outer} ({inner})")
}

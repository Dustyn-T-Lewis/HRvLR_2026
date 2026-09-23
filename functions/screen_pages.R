# The two screen pages every feature level's packet carries, drawn the same
# way so proteins, pathways and modules read side by side.

pacman::p_load(here, dplyr, ggplot2)

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
        "pROC AUC and Wilcoxon p, paired within arm; %d features",
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

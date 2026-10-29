library(here)
library(dplyr)
library(ggplot2)
source(here("04_Figures", "style", "figure_style.R"))

inputs <- c(variance = "02_HRvLR/04_Variance/c_data/variance.xlsx")
diffvar <- read_sheet(inputs[["variance"]], "diffvar")
biopsies <- read_sheet(inputs[["variance"]], "pc_scores")
per_time <- count(biopsies, timepoint, name = "n_biopsies")
per_group <- n_distinct(biopsies$subject) / n_distinct(biopsies$group)

calls <- c("HR more variable", "LR more variable", "not significant")
axis_limit <- 4
panel_data <- diffvar |>
  mutate(
    log_var_ratio = as.numeric(log_var_ratio),
    call = factor(case_when(
      p >= alpha_nominal ~ calls[3], log_var_ratio < 0 ~ calls[1],
      .default = calls[2]
    ), calls),
    beyond_axis = abs(log_var_ratio) > axis_limit
  )
counts <- panel_data |>
  summarise(
    n_tested = n(),
    n_hr = sum(call == calls[1]), n_lr = sum(call == calls[2]),
    n_fdr = sum(fdr < alpha_fdr), n_beyond = sum(beyond_axis),
    .by = timepoint
  ) |>
  mutate(
    text = sprintf(
      "%d proteins\np < %s: HR %d, LR %d\nchance %.0f each; FDR < %s: %d\nbeyond axis: %d",
      n_tested, alpha_nominal, n_hr, n_lr, alpha_nominal * n_tested / 2, alpha_fdr, n_fdr, n_beyond
    )
  )
stopifnot(
  all(is.finite(panel_data$log_var_ratio)), setequal(counts$timepoint, per_time$timepoint),
  n_distinct(per_time$n_biopsies) == 1,
  !anyDuplicated(select(panel_data, timepoint, uniprot_id))
)

panel <- panel_data |>
  filter(!beyond_axis) |>
  ggplot(aes(log_var_ratio, fill = call)) +
  geom_histogram(binwidth = 0.2, boundary = 0, linewidth = 0) +
  geom_vline(xintercept = 0, linewidth = 0.25) +
  geom_label(
    aes(-Inf, Inf, label = text), counts,
    hjust = 0, vjust = 1, size = 1.85, lineheight = 0.95, colour = "grey20",
    fill = "white", linewidth = 0, label.padding = unit(0.8, "mm"), inherit.aes = FALSE
  ) +
  scale_fill_manual(
    values = c(group_colours[["HR"]], group_colours[["LR"]], "grey88") |> rlang::set_names(calls),
    name = NULL
  ) +
  scale_x_continuous(limits = c(-axis_limit, axis_limit)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.5))) +
  facet_wrap(~timepoint, nrow = 1) +
  labs(
    x = "log variance ratio, LR relative to HR (> 0: LR people spread more)", y = "proteins",
    title = "Does one group vary more? (exploratory)",
    subtitle = sprintf(
      "DiffVar per time point, proteins observed in all %d biopsies, %d people per group",
      per_time$n_biopsies[1], per_group
    )
  ) +
  figure_theme +
  theme(
    legend.position = "bottom", strip.text = element_text(face = "bold", size = 6.5),
    panel.grid.major.x = element_blank()
  )

save_panel(panel, "S1", "c_diffvar", 125, 72)

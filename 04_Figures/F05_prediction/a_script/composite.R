# Assemble the F05 prediction composite from the leaf workbooks. BH is applied
# once here, across each figure grid (classification and continuous separately),
# because the correction spans leaves. The composite states the honest null on
# its face: AUC near chance and Q^2 <= 0 with non-significant permutation p.

pacman::p_load(here, dplyr, purrr, tidyr, ggplot2, openxlsx, patchwork)
source(here("functions", "shared_style.R"))
source(here("functions", "shared_prediction.R"))
source(here("functions", "pred_panels.R"))

root <- here("04_Figures", "F05_prediction")

read_metrics <- function(dir) {
  wb <- file.path(dir, "c_data", "source_data.xlsx")
  if (!file.exists(wb)) {
    return(NULL)
  }
  read.xlsx(wb, sheet = "metrics")
}

class_dirs <- list.dirs(file.path(root, "group_HRvLR"), recursive = TRUE) |>
  keep(~ file.exists(file.path(.x, "c_data", "source_data.xlsx")))
cont_dirs <- list.dirs(file.path(root, "phenotype"), recursive = TRUE) |>
  keep(~ file.exists(file.path(.x, "c_data", "source_data.xlsx")))

class_metrics <- map_dfr(class_dirs, read_metrics) |>
  mutate(
    q_bh = p.adjust(perm_p, "BH"),
    contrast = factor(
      contrast,
      levels = c("T1", "T2", "T3", "training", "acute")
    ),
    space = space_factor(space),
    model = factor(model, levels = c("glmnet", "spls", "pam")),
    row = interaction(model, space, sep = " / "),
    lab = sprintf("%.2f\n%s", estimate, ifelse(q_bh < 0.05, "*", ""))
  )

cont_metrics <- map_dfr(cont_dirs, read_metrics) |>
  mutate(
    q_bh = p.adjust(perm_p_q2, "BH"),
    space = space_factor(space),
    model = factor(model, levels = c("glmnet", "spls")),
    row = interaction(model, space, sep = " / "),
    lab = sprintf("%.2f\n%s", q2, ifelse(q_bh < 0.05, "*", ""))
  )

p_class <- ggplot(class_metrics, aes(contrast, row, fill = estimate)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = lab), size = 2.3, lineheight = 0.8) +
  scale_fill_gradient2(
    midpoint = 0.5, low = "#B2182B", mid = "grey92", high = "#2166AC",
    limits = c(0, 1), name = "LOSO AUC"
  ) +
  labs(
    x = NULL, y = NULL, title = "Classification (HR vs LR)",
    subtitle = paste(
      "T1 = forecast; T2/T3 = separability; training/acute = trajectory.",
      "* BH q<.05 on the permutation p."
    )
  ) +
  FIG_THEME +
  theme(
    plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE),
    axis.text.y = element_text(size = 7)
  )

cont_metrics <- cont_metrics |> mutate(q2_disp = pmin(pmax(q2, -0.6), 0.6))
p_cont <- ggplot(cont_metrics, aes(outcome, row, fill = q2_disp)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = lab), size = 2.2, lineheight = 0.8) +
  scale_fill_gradient2(
    midpoint = 0, low = "#B2182B", mid = "grey92", high = "#2166AC",
    name = expression(Q^2)
  ) +
  labs(
    x = NULL, y = NULL, title = "Continuous adaptation (baseline features)",
    subtitle = "Q^2 <= 0 means no better than the outcome mean. * BH q<.05."
  ) +
  FIG_THEME +
  theme(
    plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE),
    axis.text.x = element_text(angle = 30, hjust = 1, size = 7),
    axis.text.y = element_text(size = 7)
  )

composite <- (p_class / p_cont) +
  plot_layout(heights = c(1, 0.8)) +
  plot_annotation(
    title = "F05 — out-of-sample prediction of training response",
    subtitle = sprintf(
      paste(
        "Nested leave-one-subject-out, train-only scaling, %d-permutation",
        "null. Penalized/sparse models only. singscore is leakage-free;",
        "proteins (cohort-imputed) and modules (cohort-relative) optimistic."
      ),
      N_PERM
    ),
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(face = "bold", size = FIG_TITLE_SIZE + 1),
      plot.subtitle = element_text(face = "italic", colour = "grey30")
    )
  )

dir.create(file.path(root, "b_reports"), showWarnings = FALSE)
dir.create(file.path(root, "c_data"), showWarnings = FALSE)
save_panel(composite, file.path(root, "b_reports", "F05_composite"),
  width = 240, height = 250
)

wb <- createWorkbook()
addWorksheet(wb, "classification")
writeData(wb, "classification", class_metrics |> dplyr::select(-row, -lab))
addWorksheet(wb, "continuous")
writeData(wb, "continuous", cont_metrics |> dplyr::select(-row, -lab))
saveWorkbook(wb, file.path(root, "c_data", "F05_metrics.xlsx"),
  overwrite = TRUE
)

message(
  "F05 classification cells clearing BH q<.05: ",
  sum(class_metrics$q_bh < 0.05)
)
message(
  "F05 continuous cells clearing BH q<.05: ",
  sum(cont_metrics$q_bh < 0.05)
)

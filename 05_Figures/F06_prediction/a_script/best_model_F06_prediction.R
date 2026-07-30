suppressWarnings(suppressMessages({
  library(here)
  source(here("functions", "sweep_best_model.R"))
}))

won <- lapply(
  c("modules", "pathways", "proteins"),
  function(lv) build_best_model_level("F06_prediction", lv)
) |>
  dplyr::bind_rows()

print(as.data.frame(won[, c(
  "level", "config", "outcome", "model", "q2", "perm_p_q2", "q", "rule",
  "n_screened"
)]))

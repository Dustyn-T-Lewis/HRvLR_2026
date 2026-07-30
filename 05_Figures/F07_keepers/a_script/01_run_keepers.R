# F07: the three results figures on one sheet.
#
# A is F04 - how close each of the nine contrasts came to a BH survivor, at all
# three feature levels. B and C are F05 and F06 - the best cell every model
# reached, against that cell's own permutation null and the number of cells the
# winner was drawn from. Nothing is refitted here; all three read the workbooks
# their own stages wrote.

pacman::p_load(here, patchwork, ggplot2, openxlsx, dplyr)
source(here("functions", "keeper_panel.R"))

stem <- here("05_Figures", "F07_keepers")
for (sub in c("b_reports", "c_data")) {
  dir.create(file.path(stem, sub), recursive = TRUE, showWarnings = FALSE)
}

assoc <- association_panel()
panels <- lapply(KEEPER_ROOTS, keeper_panel)
rows <- lapply(KEEPER_ROOTS, keeper_rows)
names(rows) <- vapply(KEEPER_ROOTS, `[[`, character(1), "root")
rows$association <- association_rows()

cat(sprintf(
  "\nAssociation: %d of %d contrast x level cells clear BH .05\n",
  sum(rows$association$n_bh > 0), nrow(rows$association)
))

for (i in seq_along(KEEPER_ROOTS)) {
  r <- rows[[i]]
  cat(sprintf(
    "\n%s: %d rows, %d clear null and baseline\n",
    KEEPER_ROOTS[[i]]$label, nrow(r), sum(r$lead, na.rm = TRUE)
  ))
  print(as.data.frame(
    r |>
      arrange(desc(.data$metric)) |>
      transmute(
        .data$level, .data$model, .data$config, .data$n_screened,
        metric = round(.data$metric, 3), perm_p = signif(.data$perm_p, 3),
        null_hi = round(.data$null_hi, 3), .data$lead
      ) |>
      head(6)
  ))
}

composite <- assoc / panels[[1]] / panels[[2]] +
  plot_layout(heights = c(0.75, 1, 1)) +
  plot_annotation(tag_levels = "A")

save_panel(
  composite, file.path(stem, "b_reports", "MAIN_F07_keepers"),
  width = 190, height = 300
)

write.xlsx(rows, file.path(stem, "c_data", "F07_keepers_source_data.xlsx"))

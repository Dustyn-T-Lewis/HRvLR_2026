# F03 A: sets at FDR < 0.05 per contrast, fgsea before and after collapse, and fry.
source(here::here("05_Figures", "theme.R"))

methods <- c(n_fdr_fgsea = "fgsea", n_collapsed = "fgsea, collapsed", n_fdr_fry = "fry")
data <- read_sheet(
  "03_Pathway_Enrichment/01_run_fgsea_and_fry/c_data/01_run_fgsea_and_fry.xlsx", "set_summary"
) |>
  select(contrast, all_of(names(methods))) |>
  pivot_longer(-contrast, names_to = "method", values_to = "n_sets") |>
  mutate(
    method = factor(methods[method], methods), contrast = factor(contrast, rev(contrast_order))
  )

plot <- ggplot(data, aes(n_sets, contrast, fill = method)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  scale_fill_manual(values = c("grey75", "grey35", "#D95F02"), name = NULL) +
  labs(title = "Sets at FDR < 0.05", x = "sets", y = NULL) +
  theme_figure() +
  theme(legend.position = "bottom")
save_panel(plot, "F03", "A_set_counts")

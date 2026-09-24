# F03 C: the ten strongest collapse survivors on the primary contrast.
source(here::here("05_Figures", "theme.R"))

data <- read_sheet(
  "03_Pathway_Enrichment/01_run_fgsea_and_fry/c_data/01_run_fgsea_and_fry.xlsx", "significant"
) |>
  filter(contrast == "Training_Interaction", method == "fgsea", main) |>
  slice_min(fdr, n = 10, with_ties = FALSE) |>
  mutate(label = enrichVolcano::ev_clean_label(pathway) |> gsub(pattern = "\n", replacement = " "))

plot <- ggplot(data, aes(nes, reorder(label, nes), size = n_proteins, colour = collection)) +
  geom_vline(xintercept = 0, colour = "grey75") +
  geom_point() +
  scale_size_continuous(range = c(1.2, 3.5), name = "proteins") +
  scale_colour_brewer(palette = "Dark2", name = NULL) +
  labs(title = "Training_Interaction, fgsea", x = "NES", y = NULL) +
  theme_figure()
save_panel(plot, "F03", "C_top_sets", width = 178, height = 60)

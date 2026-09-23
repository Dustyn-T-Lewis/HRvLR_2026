# The protein packet: contrasts, their calibration, then the two screens.
# Reads only 02_Proteins c_data. Pages describe; none of them concludes.

pacman::p_load(here, dplyr, tidyr, readr, purrr, ggplot2)

source(here("functions", "contrasts.R"))
source(here("functions", "screen_pages.R"))

DEP_DIR <- here("02_Proteins", "01_Differential", "c_data")
SCR_DIR <- here("02_Proteins", "02_Classify_Associate", "c_data")

dep <- file.path(DEP_DIR, "01_dep_results.csv") |>
  read_csv(show_col_types = FALSE) |>
  mutate(contrast = factor(.data$contrast, levels = CONTRAST_NAMES))
screens <- readRDS(file.path(SCR_DIR, "protein_screens.rds"))

counts <- dep |>
  summarise(
    tested = sum(!is.na(.data$P.Value)),
    nominal = sum(.data$P.Value < 0.05, na.rm = TRUE),
    pi_score = sum(.data$sig_pi != 0L, na.rm = TRUE),
    BH = sum(.data$adj.P.Val < 0.05, na.rm = TRUE),
    .by = "contrast"
  ) |>
  mutate(expected = 0.05 * .data$tested)

p_counts <- counts |>
  pivot_longer(c("nominal", "pi_score", "BH"), names_to = "layer") |>
  mutate(layer = factor(.data$layer, c("nominal", "pi_score", "BH"))) |>
  ggplot(aes(.data$value, .data$contrast, fill = .data$layer)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_point(
    data = counts, aes(.data$expected, .data$contrast),
    inherit.aes = FALSE, shape = 124, size = 6
  ) +
  scale_fill_manual(values = c(
    nominal = "grey70", pi_score = "#4393C3", BH = "#B2182B"
  )) +
  scale_y_discrete(limits = rev) +
  labs(
    title = "Proteins called per contrast",
    subtitle = sprintf(
      "limma via proteoDA, ~ 0 + group + (1 | subject); %d proteins; %s",
      n_distinct(dep$uniprot_id), "BH within contrast"
    ),
    x = "Proteins", y = NULL, fill = "Layer",
    caption = caption(
      "Bars: proteins at nominal p < 0.05, at pi-score < 0.05 (p^|log2FC|, ",
      "Xiao 2014; no FDR control) and at BH < 0.05. Vertical tick: the count ",
      "expected at p < 0.05 if no protein responded (0.05 x proteins tested). ",
      "Data: 02_Proteins/01_Differential/c_data/01_dep_results.csv."
    )
  ) +
  FIG_THEME

p_volcano <- dep |>
  filter(!is.na(.data$P.Value)) |>
  mutate(direction = case_when(
    .data$sig_pi == 1L ~ "Up", .data$sig_pi == -1L ~ "Down", TRUE ~ "NS"
  )) |>
  ggplot(aes(.data$logFC, -log10(.data$P.Value), colour = .data$direction)) +
  geom_point(size = 0.5, alpha = 0.7) +
  facet_wrap(~contrast, nrow = 3) +
  scale_colour_manual(values = DIR_COLORS) +
  labs(
    title = "Volcano per contrast",
    subtitle = "log2 fold change against nominal p, coloured by pi-score call",
    x = "log2 fold change", y = "-log10 p", colour = "Pi-score",
    caption = caption(
      "Each point is one protein in one contrast. Colour: up or down at ",
      "pi-score < 0.05, grey otherwise. Proteins at BH < 0.05 across all ",
      "contrasts: ", sum(counts$BH), ". ",
      "Data: 02_Proteins/01_Differential/c_data/01_dep_results.csv."
    )
  ) +
  FIG_THEME

p_hist <- dep |>
  filter(!is.na(.data$P.Value)) |>
  ggplot(aes(.data$P.Value)) +
  geom_histogram(breaks = seq(0, 1, 0.05), fill = "grey55", colour = "white") +
  geom_hline(
    data = counts, aes(yintercept = .data$tested / 20),
    linetype = "dashed", colour = "#B2182B"
  ) +
  facet_wrap(~contrast, nrow = 3) +
  labs(
    title = "P-value distribution per contrast",
    subtitle = "20 bins of width 0.05",
    x = "Nominal p", y = "Proteins",
    caption = caption(
      "A uniform histogram is what a contrast with no signal produces; a ",
      "spike near zero is signal, a hump near one is over-dispersion or a ",
      "misspecified variance. Dashed line: the uniform height (tested / 20). ",
      "Data: 02_Proteins/01_Differential/c_data/01_dep_results.csv."
    )
  ) +
  FIG_THEME

SCREENS <- "02_Proteins/02_Classify_Associate/c_data/protein_screens.xlsx"
p_auc <- auc_page(screens$classify, "Protein", SCREENS)
p_assoc <- association_page(screens$chance_associate, "Protein", SCREENS)

# Gene symbols label the hit pages; a symbol two proteins share keeps its
# accession so the rows stay distinct.
gene_label <- dep |>
  distinct(.data$uniprot_id, .data$gene) |>
  mutate(label = if_else(
    is.na(.data$gene) | duplicated(.data$gene) |
      duplicated(.data$gene, fromLast = TRUE),
    paste0(coalesce(.data$gene, "?"), " (", .data$uniprot_id, ")"),
    .data$gene
  ))
label_of <- \(id) gene_label$label[match(id, gene_label$uniprot_id)]

contrast_hits <- hit_pages(
  dep |>
    transmute(
      label = label_of(.data$uniprot_id), column = .data$contrast,
      effect = .data$logFC, p = .data$P.Value, bh = .data$adj.P.Val
    ),
  CONTRAST_NAMES,
  title = "Protein contrast hits",
  subtitle = "limma via proteoDA, nominal p per contrast",
  effect_label = "log2 FC",
  data_note = "02_Proteins/01_Differential/c_data/01_dep_results.xlsx"
)

write_packet(
  c(
    list(
      "Proteins called per contrast" = p_counts,
      "Volcano per contrast" = p_volcano,
      "P-value distribution per contrast" = p_hist
    ),
    list_flatten(
      list("Protein contrast hits" = contrast_hits),
      name_spec = "{outer} ({inner})"
    ),
    list(
      "Protein AUC per classification task" = p_auc,
      "Protein association with phenotype" = p_assoc
    ),
    screen_hit_pages(screens, "Protein", SCREENS, label_of)
  ),
  here("02_Proteins", "03_Packet", "b_reports", "02_Proteins_packet.pdf"),
  title = "HRvLR 02 Proteins"
)

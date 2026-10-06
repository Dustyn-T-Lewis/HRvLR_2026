# Figure 1: the design and the whole-proteome picture before any single result. Reads saved
# outputs only. The one computation is the phenotype PCA, a projection with no test; every
# p-value plotted comes from 02_Contrasts or 03_Contrasts.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})
source(here("R", "helpers.R"))
source(here("R", "figure_theme.R"))
source(here("R", "figure_panels.R"))

stage <- here("02_HRvLR_groups", "Figures", "F01_Overview")
figure_dir <- file.path(stage, "b_reports")
out <- file.path(stage, "c_data")
walk(c(figure_dir, out), dir.create, recursive = TRUE, showWarnings = FALSE)

inputs <- c(
  contrasts = "02_HRvLR_groups/02_DEP/c_data/02_contrasts.xlsx",
  set_tests = "02_HRvLR_groups/03_Pathways/c_data/set_tests.rds",
  gene_sets = "01_Trunk/04_Gene_Sets/c_data/gene_sets.rds",
  phenotype = "00_Input/phenotype.csv",
  design = "02_HRvLR_groups/Figures/F01_Overview/assets/design.png"
)
paths <- map_chr(inputs, here)
if (!all(file.exists(paths))) {
  stop("Missing inputs: ", paste(inputs[!file.exists(paths)], collapse = ", "))
}
read_sheet <- \(sheet) readxl::read_excel(paths[["contrasts"]], sheet)
set_tests <- readRDS(paths[["set_tests"]])
proteins <- set_tests$protein_results
gene_sets <- readRDS(paths[["gene_sets"]])

# Panel a: the BioRender design panel, pinned under its title so it sits level with b and c.
panel_design <- ggplot() +
  annotation_custom(grid::rasterGrob(
    png::readPNG(paths[["design"]]),
    y = grid::unit(1, "npc"), just = "top", interpolate = TRUE
  )) +
  labs(
    title = "Design and contrasts",
    subtitle = "16 adults, 45 biopsies; contrast colours carry through b to f"
  ) +
  theme_void(base_size = 7) +
  theme(
    plot.title = figure_theme$plot.title, plot.subtitle = figure_theme$plot.subtitle,
    plot.title.position = "plot"
  )

# Panel b: proteins per contrast. The light bar counts every protein at p < 0.05 (or Pi < 0.05);
# the dark bar, drawn over it from zero, counts the Pi < 0.05 subset. Each row's name sits in a
# box in its family's tint and border, as in panel a.
contrast_labels <- tribble(
  ~contrast, ~group, ~label,
  "Training_HR", "Training", "HR",
  "Training_LR", "Training", "LR",
  "Acute_HR", "Acute", "HR",
  "Acute_LR", "Acute", "LR",
  "Baseline_HRvLR", "HR - LR", "T1",
  "Trained_HRvLR", "HR - LR", "T2",
  "Acute_HRvLR", "HR - LR", "T3",
  "Training_Interaction", "Interaction", "Training",
  "Acute_Interaction", "Interaction", "Acute"
) |>
  mutate(group = factor(group, unique(group)), contrast = factor(contrast, rev(contrast)))
protein_counts <- proteins |>
  mutate(direction = if_else(logFC > 0, "Up", "Down")) |>
  summarise(
    pi = sum(pi_score < 0.05), nominal = sum(p < 0.05 | pi_score < 0.05),
    fdr = sum(fdr < 0.05), .by = c(contrast, direction)
  ) |>
  complete(contrast, direction, fill = list(pi = 0, nominal = 0, fdr = 0)) |>
  inner_join(mutate(contrast_labels, contrast = as.character(contrast)), by = "contrast")
pi_text <- "\u03a0 < 0.05"
bars <- protein_counts |>
  pivot_longer(c(nominal, pi), names_to = "tier", values_to = "n") |>
  mutate(
    signed = if_else(direction == "Up", n, -n),
    contrast = factor(contrast, levels(contrast_labels$contrast)),
    fill = paste(direction, if_else(tier == "pi", pi_text, "p < 0.05"))
  ) |>
  # The Pi bar is drawn last, over the nominal bar.
  arrange(tier == "pi")
bar_limit <- max(protein_counts$nominal)
bar_fills <- set_names(
  c(direction_colours[["Down"]], "#F1B6DA", "#B8E186", direction_colours[["Up"]]),
  paste(c("Down", "Down", "Up", "Up"), c(pi_text, "p < 0.05", "p < 0.05", pi_text))
)
row_boxes <- mutate(contrast_labels, family = as.character(group), x = -bar_limit * 1.06)
panel_proteins <- ggplot(bars, aes(signed, contrast)) +
  geom_col(aes(fill = fill), width = 0.7, position = "identity") +
  geom_vline(xintercept = 0, linewidth = 0.25) +
  # Families share the fill scale with the bars; the legend lists the bars only.
  geom_label(
    aes(x, contrast, label = label, fill = family, colour = family), row_boxes,
    hjust = 1, size = 2.1, fontface = "bold", text.colour = "black",
    linewidth = 0.35, label.r = unit(0.6, "mm"), label.padding = unit(0.7, "mm")
  ) +
  scale_fill_manual(
    values = c(bar_fills, contrast_tone("fill")), breaks = names(bar_fills), name = NULL
  ) +
  scale_colour_manual(values = contrast_tone("line"), guide = "none") +
  scale_x_continuous(
    labels = abs, limits = c(-bar_limit * 1.45, bar_limit), breaks = seq(-100, 100, 50)
  ) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  coord_cartesian(clip = "off") +
  labs(
    x = "proteins (down | up)", y = NULL, title = "Proteins per contrast",
    subtitle = "1,900 tested; none at FDR < 0.05"
  ) +
  figure_theme +
  theme(
    strip.text = element_blank(), axis.text.y = element_blank(),
    legend.position = "bottom", panel.grid.major.y = element_blank(),
    panel.spacing.y = unit(1.2, "mm")
  ) +
  guides(fill = guide_legend(nrow = 2))

# Panel c: proteome and phenotype PCA, each with its PERMANOVA (02_Contrasts). The phenotype is
# the six size and strength traits at T1 and T2, scaled; the one biopsy missing a leg-extension
# value is left out rather than filled. Arm is not tested on the phenotype, because HR and LR were
# split on how much these traits changed.
permanova <- read_sheet("permanova")
test_row <- \(data_, term_, window_) {
  filter(permanova, data == data_, term == term_, window == window_, arms == "both")
}
permanova_text <- \(label, test) sprintf("%s R\u00b2 %.2f, p = %.2g", label, test$r2, test$p)
time_test <- test_row("proteome", "time", "all")
arm_test <- test_row("proteome", "arm", "subject mean")
trait_time_test <- test_row("phenotype", "time", "all")
pca_proteome <- read_sheet("pca")
panel_proteome <- pca_panel(
  pca_proteome, read_sheet("pca_variance")$variance_explained,
  "Proteome", "45 biopsies, 1,900 proteins",
  paste(
    "PERMANOVA", permanova_text("time", time_test), permanova_text("arm", arm_test),
    sep = "\n"
  )
)
phenotype <- phenotype_tables(readr::read_csv(paths[["phenotype"]], show_col_types = FALSE))
traits <- filter(phenotype$biopsy, if_all(all_of(phenotype$traits), \(x) !is.na(x)))
trait_pca <- prcomp(scale(traits[phenotype$traits]))
pca_phenotype <- bind_cols(
  select(traits, subject, arm, timepoint), as_tibble(trait_pca$x[, 1:2])
)
panel_phenotype <- pca_panel(
  pca_phenotype, summary(trait_pca)$importance[2, ],
  "Phenotype", sprintf("%d biopsies, six size and strength traits", nrow(traits)),
  paste("PERMANOVA", permanova_text("time", trait_time_test), sep = "\n")
) +
  # No T3 biopsy has traits, so its key would differ from the proteome's and print twice.
  guides(colour = "none", fill = "none", shape = "none")
panel_pca <- (panel_proteome | panel_phenotype) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# Panel d: pathways at fgsea FDR < 0.05 in both arms of a window (HR - LR: at least two of the
# three time points; Interaction: either interaction, strongest five), from four collections,
# each shown in all nine contrasts. A set sharing half its measured proteins with a stronger set
# already shown is dropped. Opposite signs count only among the significant columns. fgsea calls
# 64 and 206 interaction sets and fry none, so the Interaction group stays small.
windows <- list(
  Training = c(HR = "Training_HR", LR = "Training_LR"),
  Acute = c(HR = "Acute_HR", LR = "Acute_LR"),
  `HR - LR` = c(T1 = "Baseline_HRvLR", T2 = "Trained_HRvLR", T3 = "Acute_HRvLR"),
  Interaction = c(Training = "Training_Interaction", Acute = "Acute_Interaction")
)
needed <- c(Training = 2, Acute = 2, `HR - LR` = 2, Interaction = 1)
row_cap <- c(Training = Inf, Acute = Inf, `HR - LR` = Inf, Interaction = 5)
collections <- c(
  Hallmark = "#332288", Reactome = "#88CCEE", KEGG_Legacy = "#DDCC77", GOBP = "#777777"
)
members <- map(gene_sets$sets, \(s) intersect(s, gene_sets$gene_universe))
fgsea <- set_tests$set_tests |>
  filter(method == "fgsea", collection %in% names(collections)) |>
  select(contrast, set_id, collection, pathway, nes, fdr, leading_edge)
select_window <- function(name) {
  rows <- filter(fgsea, contrast %in% windows[[name]])
  eligible <- rows |>
    summarise(
      hits = sum(fdr < 0.05),
      same_sign = n_distinct(sign(nes[fdr < 0.05])) == 1,
      strength = sum(abs(nes)),
      # Several Reactome signalling sets carry only proteasome subunits; the label says so.
      proteasome = mean(startsWith(as.character(unique(unlist(leading_edge[fdr < 0.05]))), "PSM")),
      .by = c(set_id, pathway)
    ) |>
    filter(hits >= needed[[name]]) |>
    arrange(desc(same_sign), desc(strength))
  kept <- head(prune_overlapping(eligible$set_id, members), row_cap[[name]])
  eligible |>
    filter(set_id %in% kept) |>
    transmute(
      facet = name, set_id, order = match(set_id, kept),
      label = paste0(
        if_else(proteasome >= 0.5, "Proteasome: ", ""), pathway_label(pathway),
        if_else(same_sign, "", " *")
      )
    )
}
columns <- imap(windows, \(contrasts, window) {
  tibble(window, arm = names(contrasts), contrast = unname(contrasts))
}) |>
  list_rbind() |>
  mutate(
    column = factor(paste(window, arm), c(" ", paste(window, arm))),
    markup = in_colour(arm, arm_shade(window, arm)),
    window = factor(window, c("", names(windows)))
  )
pathways <- map(names(windows), select_window) |>
  list_rbind() |>
  inner_join(fgsea, by = "set_id", relationship = "many-to-many") |>
  inner_join(columns, by = "contrast") |>
  mutate(
    facet = factor(facet, names(windows)), row = paste(facet, set_id),
    log_fdr = pmin(-log10(fdr), 12)
  )
# A narrow first column marks each set's database.
databases <- distinct(pathways, row, facet, order, collection) |>
  mutate(window = factor("", levels(columns$window)), column = factor(" ", levels(columns$column)))
nes_limit <- ceiling(max(abs(pathways$nes)) * 10) / 10
panel_pathways <- pathways |>
  mutate(row = reorder(row, -order)) |>
  ggplot(aes(column, row)) +
  geom_point(aes(size = log_fdr, fill = nes), shape = 21, stroke = 0.2) +
  geom_point(
    aes(colour = collection), mutate(databases, row = reorder(row, -order)),
    shape = 15, size = 1.8
  ) +
  scale_colour_manual(
    values = collections, breaks = names(collections),
    labels = c("Hallmark", "Reactome", "KEGG", "GO BP"), name = "Database"
  ) +
  scale_y_discrete(labels = deframe(distinct(pathways, row, label))) +
  scale_x_discrete(labels = c(` ` = "", deframe(select(columns, column, markup)))) +
  scale_fill_distiller(
    palette = "PiYG", direction = 1, limits = c(-nes_limit, nes_limit), name = "NES"
  ) +
  scale_size_area(
    max_size = 2.8, limits = c(0, 12), breaks = c(2, 6, 12),
    labels = c("2", "6", "12+"), name = "-log10 FDR"
  ) +
  guides(
    fill = guide_colourbar(order = 1), size = guide_legend(order = 2),
    colour = guide_legend(order = 3)
  ) +
  ggh4x::facet_grid2(
    facet ~ window,
    scales = "free", space = "free", switch = "y",
    strip = contrast_strips(x = levels(columns$window), y = names(windows))
  ) +
  labs(
    x = NULL, y = NULL, title = "Pathways shared within each window",
    subtitle = paste(
      "fgsea FDR < 0.05 in both arms of a window (interaction: either, top 5),",
      "shown in all nine contrasts; * opposite signs"
    )
  ) +
  figure_theme +
  theme(
    axis.text.y = element_text(size = 5.2),
    axis.text.x.bottom = ggtext::element_markdown(),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 90, size = 5.5, face = "bold"),
    strip.text.x = element_text(size = 6, face = "bold"),
    panel.spacing = unit(1, "mm")
  )

# Panels e and f: each arm's p < 0.05 sets for a window placed in the other arm's ranking. The Pi
# sets and the training-to-bout trajectories are in S03_Signatures.
fry_rows <- read_sheet("fry_concordance")
barcodes <- map(
  set_names(c("Training", "Acute")),
  \(window) signature_barcode(proteins, fry_rows, window)
)

# Every plotted number traces to its source table before anything is written.
contrast_summary <- read_sheet("contrast_summary")
counted <- summarise(protein_counts, fdr = sum(fdr), .by = contrast) |>
  inner_join(contrast_summary, by = "contrast")
barcode_numbers <- map(barcodes, "labels") |> list_rbind()
fry_numbers <- fry_rows |>
  filter(startsWith(set, "nominal_"), sub("_.*", "", set_from) == sub("_.*", "", ranked_on)) |>
  transmute(
    window = sub("_.*", "", set_from), chosen = sub(".*_", "", set_from),
    ranked_in = sub(".*_", "", ranked_on), direction = sub(".*_", "", set), fdr_source = fdr
  )
traced <- inner_join(
  barcode_numbers, fry_numbers,
  by = c("window", "chosen", "ranked_in", "direction")
)
stopifnot(
  nrow(counted) == 9, all(counted$fdr == counted$n_fdr_05),
  nrow(pca_proteome) == 45, nrow(pca_phenotype) == 31,
  nrow(time_test) == 1, nrow(arm_test) == 1, nrow(trait_time_test) == 1,
  all(count(pathways, row)$n == nrow(columns)),
  nrow(traced) == 8, all(traced$fdr == traced$fdr_source)
)

rows <- list(
  # free() keeps a and c from aligning to the long labels in panel d.
  (free(panel_design) | panel_proteins | free(panel_pca)) +
    plot_layout(widths = c(1.1, 1.45, 2)),
  (panel_pathways | barcodes$Training$plot | barcodes$Acute$plot) +
    plot_layout(widths = c(1.45, 1, 1))
)
figure <- wrap_plots(rows, ncol = 1, heights = c(66, 98)) +
  # The phenotype PCA shares tag c with the proteome PCA.
  plot_annotation(tag_levels = list(c("a", "b", "c", "", "d", "e", "f"))) &
  figure_tags
save_figure(figure, figure_dir, "f01_overview", width = 300, height = 175)

sheets <- list(
  protein_counts = select(protein_counts, group, label, contrast, direction, pi, nominal, fdr),
  pca_proteome = pca_proteome,
  pca_phenotype = pca_phenotype,
  permanova = permanova,
  pathways = transmute(
    pathways,
    selected_in = as.character(facet), window = as.character(window),
    column = as.character(arm), set_id, collection, label, nes, fdr
  ),
  barcode = select(barcode_numbers, window, chosen, ranked_in, direction, n, same_sign_share, fdr)
)
holds <- c(
  "Panel b: per contrast and direction, proteins at Pi < 0.05, at p < 0.05 or Pi, at FDR < 0.05.",
  "Panel c, left: PC1 and PC2 per sample (02_Contrasts, pca).",
  "Panel c, right: PC1 and PC2 per biopsy from the six scaled size and strength traits.",
  "Panel c: PERMANOVA for the proteome (time, arm) and the phenotype (time) (02_Contrasts).",
  "Panel d: each selected set in all nine contrasts, with NES and fgsea FDR.",
  "Panels e, f: each arm's p < 0.05 sets in the other arm's ranking: n, % same way, fry FDR."
)
write_workbook(
  sheets, holds, file.path(out, "f01_overview.xlsx"), inputs, paths,
  note = "Panel a is the BioRender design panel (assets/design.png); it reads no data."
)
message("wrote f01_overview.pdf, f01_overview.png and f01_overview.xlsx")

# The pathway packet: the set universe and its themes, the set tests of the
# nine contrasts, per-sample scores, then the two screens.
# Reads only 03_Pathways c_data. Pages describe; none of them concludes.

pacman::p_load(here, dplyr, tidyr, forcats, purrr, ggplot2, patchwork)

source(here("functions", "contrasts.R"))
source(here("functions", "screen_pages.R"))

stage <- function(...) here("03_Pathways", ...)
gs <- readRDS(stage("01_Gene_Sets", "c_data", "gene_sets.rds"))
st <- readRDS(stage("02_Set_Tests", "c_data", "set_tests.rds"))
ss <- readRDS(stage("03_Set_Scores", "c_data", "set_scores.rds"))
screens <- readRDS(stage(
  "04_Classify_Associate", "c_data", "pathway_screens.rds"
))
tests <- st$set_tests

p_collections <- count(gs$catalog, .data$collection) |>
  ggplot(aes(.data$n, fct_reorder(.data$collection, .data$n))) +
  geom_col(fill = "grey55") +
  geom_text(aes(label = .data$n), hjust = -0.2, size = 3) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(x = "Sets kept", y = NULL) +
  FIG_THEME
p_themes <- gs$catalog |>
  filter(!is.na(.data$theme)) |>
  count(.data$theme) |>
  ggplot(aes(.data$n, fct_reorder(.data$theme, .data$n))) +
  geom_col(fill = DB_COLORS[["GO:BP"]]) +
  labs(x = "GO:BP sets", y = NULL) +
  FIG_THEME +
  theme(axis.text.y = element_text(size = 6.5))
p_catalog <- p_collections + p_themes + plot_layout(widths = c(1, 2)) +
  plot_annotation(
    title = "Gene-set universe and GO-Slim themes",
    subtitle = sprintf(
      "MSigDB %s and GO-Slim generic (%s); 15 to 500 detected members per set",
      gs$msigdb_release, gs$goslim_release
    ),
    caption = caption(
      "Left: sets kept per collection after the detected-member floor. Right: ",
      "GO:BP sets grouped by their most specific GO-Slim ancestor (",
      sum(!is.na(gs$catalog$theme)), " of ",
      sum(gs$catalog$collection == "GO:BP"),
      " GO:BP sets fall under a slim term; the rest carry no theme). ",
      "Data: 03_Pathways/01_Gene_Sets/c_data/gene_sets.xlsx, sheets catalog ",
      "and themes."
    ),
    theme = FIG_THEME
  )

fry_counts <- tests |>
  summarise(
    nominal = sum(.data$fry_p < 0.05),
    ratio = .data$nominal / (0.05 * n()),
    fdr = sum(.data$fry_fdr < 0.05),
    .by = c("contrast", "collection")
  )
p_fry <- ggplot(fry_counts, aes(.data$collection, .data$contrast)) +
  geom_tile(aes(fill = .data$ratio), colour = "white") +
  geom_text(aes(label = paste(.data$nominal, "/", .data$fdr)), size = 3) +
  scale_fill_gradient2(
    low = "#4393C3", mid = "white", high = "#D6604D", midpoint = 1
  ) +
  scale_y_discrete(limits = rev) +
  labs(
    title = "fry set tests per contrast and collection",
    subtitle = sprintf(
      "limma::fry, protein design, subject block, cor %.3f; BH per %s",
      st$correlation, "contrast x collection"
    ),
    x = NULL, y = NULL, fill = "Nominal /\nexpected",
    caption = caption(
      "Label: sets at nominal fry p < 0.05 / sets at BH < 0.05. Fill: nominal ",
      "sets over the count expected by chance (0.05 x sets in the ",
      "collection); white is what no signal looks like. ",
      "Data: 03_Pathways/02_Set_Tests/c_data/set_tests.xlsx, sheet summary."
    )
  ) +
  FIG_THEME

shown <- tests |>
  filter(.data$main, .data$fgsea_padj < 0.05) |>
  summarise(best = min(.data$fgsea_padj), .by = "set") |>
  slice_min(.data$best, n = 40, with_ties = FALSE) |>
  pull("set")
p_nes <- tests |>
  filter(.data$set %in% shown) |>
  mutate(
    group = coalesce(.data$theme, .data$collection),
    label = clean_set_name(.data$set)
  ) |>
  ggplot(aes(.data$contrast, .data$label, fill = .data$nes)) +
  geom_tile(colour = "white") +
  geom_point(
    data = \(d) filter(d, .data$fry_fdr < 0.05), size = 1.2, shape = 16
  ) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B") +
  labs(
    title = "Leading non-redundant sets by fgsea",
    subtitle = sprintf(
      "fgsea on moderated t; collapsePathways main sets, padj < 0.05; top %d",
      length(shown)
    ),
    x = NULL, y = NULL, fill = "NES",
    caption = caption(
      "Rows: sets fgsea calls at padj < 0.05 in at least one contrast and ",
      "that collapsePathways keeps as non-redundant, grouped by GO-Slim theme ",
      "(or collection where a set has none). Fill: NES in every contrast. ",
      "Dot: fry BH < 0.05 in that contrast. fgsea's gene-permutation p ",
      "assumes exchangeable proteins, so fry is the test to quote. ",
      "Data: 03_Pathways/02_Set_Tests/c_data/set_tests.xlsx."
    )
  ) +
  FIG_THEME +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    axis.text.y = element_text(size = 6),
    strip.text.y = element_text(angle = 0, size = 6.5, hjust = 0)
  )

theme_nes <- tests |>
  filter(!is.na(.data$theme)) |>
  summarise(
    nes = mean(.data$nes, na.rm = TRUE),
    n = n(),
    .by = c("theme", "contrast")
  )
theme_order <- theme_nes |>
  pivot_wider(
    id_cols = "theme", names_from = "contrast", values_from = "nes"
  ) |>
  tibble::column_to_rownames("theme") |>
  as.matrix() |>
  (\(m) rownames(m)[hclust(dist(m))$order])()
theme_nes$theme <- factor(theme_nes$theme, levels = theme_order)
p_theme <- theme_nes |>
  ggplot(aes(.data$contrast, .data$theme, fill = .data$nes)) +
  geom_tile(colour = "white") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B") +
  labs(
    title = "Theme-level direction per contrast",
    subtitle = "Mean fgsea NES of the GO:BP sets in each theme; rows clustered",
    x = NULL, y = NULL, fill = "Mean NES",
    caption = caption(
      "Each tile averages the NES of every GO:BP set sharing a GO-Slim ",
      "theme; it summarises direction, not significance, and a theme of ",
      "few sets moves with any one of them. Sets per theme: page 2. ",
      "Data: 03_Pathways/02_Set_Tests/c_data/set_tests.xlsx."
    )
  ) +
  FIG_THEME +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    axis.text.y = element_text(size = 6.5)
  )

pca <- prcomp(t(ss$scores))
var_pc <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
p_pca <- ss$meta |>
  mutate(pc1 = pca$x[, 1], pc2 = pca$x[, 2]) |>
  ggplot(aes(.data$pc1, .data$pc2)) +
  geom_path(aes(group = .data$subject), colour = "grey80", linewidth = 0.3) +
  geom_point(aes(colour = .data$timepoint, shape = .data$arm), size = 2.4) +
  scale_colour_manual(values = TIME_COLORS) +
  labs(
    title = "Samples in pathway-score space",
    subtitle = sprintf(
      "PCA of %d singscore set scores x %d samples",
      nrow(ss$scores), ncol(ss$scores)
    ),
    x = sprintf("PC1 (%.1f%%)", var_pc[1]),
    y = sprintf("PC2 (%.1f%%)", var_pc[2]),
    colour = "Timepoint", shape = "Arm",
    caption = caption(
      "Each point is one biopsy scored on every set; grey lines join a ",
      "subject's own biopsies. Scores are rank-based within each sample ",
      "(singscore), computed on the missForest matrix. ",
      "Data: 03_Pathways/03_Set_Scores/c_data/set_scores.csv."
    )
  ) +
  FIG_THEME

SCREENS <- "03_Pathways/04_Classify_Associate/c_data/pathway_screens.xlsx"
p_auc <- auc_page(screens$classify, "Pathway", SCREENS)
p_assoc <- association_page(
  chance_table(screens$associate, "window", "phenotype"), "Pathway", SCREENS
)

# Set labels carry a collection tag, so a Reactome and a GO:BP set with the
# same wording stay on separate rows.
COLLECTION_TAG <- c(
  Hallmark = "[H]", Reactome = "[R]", "GO:BP" = "[GO]", "GO Slim" = "[Slim]"
)
set_label <- gs$catalog |>
  transmute(
    set = .data$set,
    label = paste(
      COLLECTION_TAG[.data$collection], clean_set_name(.data$set, 55)
    )
  )
label_of <- \(set) set_label$label[match(set, set_label$set)]

fry_hits <- hit_pages(
  tests |>
    transmute(
      label = label_of(.data$set), column = .data$contrast,
      effect = .data$nes, p = .data$fry_p, bh = .data$fry_fdr
    ),
  CONTRAST_NAMES,
  title = "Pathway contrast hits (fry)",
  subtitle = "fry p per contrast, filled by fgsea NES for direction",
  effect_label = "NES",
  data_note = "03_Pathways/02_Set_Tests/c_data/set_tests.xlsx"
)

write_packet(
  c(
    list(
      "Gene-set universe and GO-Slim themes" = p_catalog,
      "fry set tests per contrast and collection" = p_fry
    ),
    list_flatten(
      list("Pathway contrast hits (fry)" = fry_hits),
      name_spec = "{outer} ({inner})"
    ),
    list(
      "Leading non-redundant sets by fgsea" = p_nes,
      "Theme-level direction per contrast" = p_theme,
      "Samples in pathway-score space" = p_pca,
      "Pathway AUC per classification task" = p_auc,
      "Pathway association with phenotype" = p_assoc
    ),
    screen_hit_pages(screens, "Pathway", SCREENS, label_of)
  ),
  stage("05_Packet", "b_reports", "03_Pathways_packet.pdf"),
  title = "HRvLR 03 Pathways"
)

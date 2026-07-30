# F02 Panel F: GO GSEA per contrast, deduplicated, split up / down (BP / CC / MF)
# GSEA on the per-contrast moderated-t ranks against the three GO ontologies. Raw
# significant-term counts are inflated by GO redundancy, so each ontology's
# significant sets (BH < 0.05) are collapsed to non-redundant representatives
# (fgsea::collapsePathways) before counting. Horizontal diverging layout: up-regulated
# right, down-regulated left, per-contrast background bands.

pacman::p_load(here, dplyr, tidyr, tibble, ggplot2, fgsea, msigdbr)

if (!exists("meta")) source(here("05_Figures", "F02_proteome", "a_script", "setup.R"))

lighten <- function(col, amount = 0.45) {
  v <- grDevices::col2rgb(col) / 255
  m <- v + (1 - v) * amount
  grDevices::rgb(m[1], m[2], m[3])
}
ONT_COLORS_LIGHT <- vapply(ONT_COLORS, lighten, character(1))

go_sets <- function(sub) {
  g <- tryCatch(
    msigdbr(species = "Homo sapiens", collection = "C5", subcollection = paste0("GO:", sub)),
    error = function(e) msigdbr(species = "Homo sapiens", category = "C5", subcategory = paste0("GO:", sub))
  )
  split(g$gene_symbol, g$gs_name)
}
pathways <- lapply(c(BP = "BP", CC = "CC", MF = "MF"), go_sets)

display_contrasts <- setdiff(MAIN_CONTRASTS, c("Training_Interaction", "Acute_Interaction"))

collapsed_counts <- function(res, sets, ranks) {
  sig <- res[!is.na(res$padj) & res$padj < 0.05, ]
  if (nrow(sig) > 1) {
    keep <- collapsePathways(sig[order(sig$pval), ], sets, ranks)$mainPathways
    sig <- sig[sig$pathway %in% keep, ]
  }
  c(Up = sum(sig$NES > 0), Down = sum(sig$NES < 0))
}

set.seed(42)
go_counts <- lapply(display_contrasts, function(cn) {
  ranks <- setNames(dep_df[[paste0("t_", cn)]], dep_df$gene)
  ranks <- ranks[!is.na(ranks) & !is.na(names(ranks)) & !duplicated(names(ranks))]
  ranks <- sort(ranks, decreasing = TRUE)
  lapply(names(ONT_COLORS), function(ont) {
    n <- collapsed_counts(fgsea(pathways[[ont]], ranks, minSize = 10, maxSize = 500), pathways[[ont]], ranks)
    tibble(contrast = cn, ontology = ont, Up = n[["Up"]], Down = n[["Down"]])
  }) |> bind_rows()
}) |>
  bind_rows()

F02_AUDIT[["panel_F_go_counts"]] <- go_counts

plot_df <- go_counts |>
  pivot_longer(c(Up, Down), names_to = "direction", values_to = "n_terms") |>
  mutate(
    signed = if_else(direction == "Down", -n_terms, n_terms),
    contrast = factor(contrast, levels = rev(display_contrasts)),
    ontology = factor(ontology, levels = names(ONT_COLORS)),
    fill_key = paste(ontology, direction, sep = "_")
  )

FILL_VALUES <- c(
  setNames(unname(ONT_COLORS), paste0(names(ONT_COLORS), "_Up")),
  setNames(unname(ONT_COLORS_LIGHT), paste0(names(ONT_COLORS), "_Down"))
)

totals <- plot_df |>
  group_by(contrast, direction) |>
  summarise(tot = sum(n_terms), .groups = "drop") |>
  mutate(signed_tot = if_else(direction == "Down", -tot, tot)) |>
  filter(tot > 0)

PH_W <- 120
PH_H <- 110

band_levels <- levels(plot_df$contrast)
band_df <- tibble(
  xmin = seq_along(band_levels) - 0.5, xmax = seq_along(band_levels) + 0.5,
  band = CONTRAST_COLORS[band_levels]
)

pF <- ggplot(plot_df, aes(contrast, signed, fill = fill_key)) +
  geom_rect(
    data = band_df, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE, fill = scales::alpha(band_df$band, 0.14),
    color = "grey75", linewidth = 0.2
  ) +
  geom_col(width = 0.7, color = "grey30", linewidth = 0.25) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "grey30") +
  geom_text(
    data = totals, aes(contrast, signed_tot, label = tot),
    inherit.aes = FALSE, hjust = ifelse(totals$direction == "Down", 1.3, -0.3),
    size = FIG_GEOM_TEXT, fontface = "bold", color = "grey15"
  ) +
  scale_fill_manual(
    values = FILL_VALUES,
    breaks = paste0(names(ONT_COLORS), "_Up"), labels = names(ONT_COLORS), name = NULL
  ) +
  scale_x_discrete(labels = CTR_SHORT) +
  scale_y_continuous(labels = abs, expand = expansion(mult = c(0.12, 0.12))) +
  coord_flip(clip = "off") +
  labs(x = NULL, y = "Down  <-  GO terms  ->  Up") +
  FIG_THEME +
  theme(
    legend.position = "bottom",
    axis.text.y = element_text(face = "bold"),
    plot.margin = margin(6, 6, 4, 4)
  )

save_png(pF, file.path(RPT_DIR, "panels", "panel_f_pathways"), PH_W, PH_H)
cat("F02 Panel F done.\n")

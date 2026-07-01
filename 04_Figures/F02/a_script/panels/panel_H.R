# F02 Panel H: GO enrichment per contrast, split up / down (BP / CC / MF, fgsea)
# GSEA on the per-contrast moderated-t ranks against the three GO ontologies.
# Significant sets (BH < 0.05) split by NES sign: up-regulated above zero (solid,
# dark), down-regulated below zero (one shade lighter, dotted outline). Stacked by
# ontology, per-contrast background bands carry the contrast code.

pacman::p_load(here, dplyr, tidyr, tibble, ggplot2, fgsea, msigdbr)

if (!exists("meta")) source(here("04_Figures", "F02", "a_script", "HRvLR_F02_setup.R"))

ONT_COLORS <- c(BP = "#66C2A5", CC = "#FC8D62", MF = "#8DA0CB")

lighten <- function(col, amount = 0.45) {
  rgb_vals <- grDevices::col2rgb(col) / 255
  mixed <- rgb_vals + (1 - rgb_vals) * amount
  grDevices::rgb(mixed[1], mixed[2], mixed[3])
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

set.seed(42)
go_counts <- lapply(display_contrasts, function(cn) {
  ranks <- setNames(dep_df[[paste0("t_", cn)]], dep_df$gene)
  ranks <- ranks[!is.na(ranks) & !is.na(names(ranks)) & !duplicated(names(ranks))]
  ranks <- sort(ranks, decreasing = TRUE)
  lapply(names(ONT_COLORS), function(ont) {
    res <- fgsea(pathways[[ont]], ranks, minSize = 10, maxSize = 500)
    sig <- res[!is.na(res$padj) & res$padj < 0.05, ]
    n_up <- sum(sig$NES > 0)
    n_dn <- sum(sig$NES < 0)
    tibble(contrast = cn, ontology = ont, Up = n_up, Down = n_dn)
  }) |> bind_rows()
}) |>
  bind_rows()

write.csv(go_counts, file.path(DAT_DIR, "audit_panel_H_go_counts.csv"), row.names = FALSE)

plot_df <- go_counts |>
  pivot_longer(c(Up, Down), names_to = "direction", values_to = "n_terms") |>
  mutate(
    signed = if_else(direction == "Down", -n_terms, n_terms),
    contrast = factor(contrast, levels = display_contrasts),
    ontology = factor(ontology, levels = names(ONT_COLORS)),
    direction = factor(direction, levels = c("Up", "Down")),
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

PH_W <- 150
PH_H <- PANEL_H

band_df <- tibble(
  contrast = factor(display_contrasts, levels = display_contrasts),
  xmin = seq_along(display_contrasts) - 0.5,
  xmax = seq_along(display_contrasts) + 0.5,
  band = CONTRAST_COLORS[display_contrasts]
)

pH <- ggplot(plot_df, aes(contrast, signed, fill = fill_key, linetype = direction)) +
  geom_rect(
    data = band_df,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE, fill = scales::alpha(band_df$band, 0.14),
    color = "grey75", linewidth = 0.2
  ) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_col(width = 0.74, color = "grey20", linewidth = 0.3) +
  geom_text(
    data = totals, aes(contrast, signed_tot, label = tot),
    inherit.aes = FALSE, vjust = ifelse(totals$direction == "Up", -0.4, 1.3),
    size = FIG_GEOM_TEXT, fontface = "bold"
  ) +
  scale_fill_manual(
    values = FILL_VALUES,
    breaks = paste0(names(ONT_COLORS), "_Up"), labels = names(ONT_COLORS),
    name = NULL
  ) +
  scale_linetype_manual(
    values = c(Up = "solid", Down = "dotted"),
    guide = guide_legend(override.aes = list(fill = "grey85"))
  ) +
  scale_x_discrete(labels = CTR_SHORT) +
  scale_y_continuous(expand = expansion(mult = c(0.10, 0.10))) +
  labs(
    title = "Pathway Enrichment (GO GSEA, Up / Down)",
    subtitle = "Significant gene sets (BH < 0.05); up above, down below zero",
    x = NULL, y = "Down  <-  GO terms  ->  Up", fill = NULL, linetype = NULL, tag = "H"
  ) +
  FIG_THEME +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 30, hjust = 1))

save_png(pH, file.path(RPT_DIR, "panels", "panel_h_pathways"), PH_W, PH_H)
cat("F02 Panel H done.\n")

# F02 Panel H: GO enrichment per contrast (BP / CC / MF, fgsea) ----
# GSEA on the per-contrast moderated-t ranks against the three GO ontologies;
# counts significant gene sets (BH < 0.05). Mirrors YvO's pathway-enrichment bar.

setwd(rprojroot::find_rstudio_root_file())
if (!exists("meta")) source("04_Figures/F02/a_script/HRvLR_F02_setup.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(fgsea)
  library(msigdbr)
})

ONT_COLORS <- c(BP = "#66C2A5", CC = "#FC8D62", MF = "#8DA0CB")

go_sets <- function(sub) {
  g <- tryCatch(
    msigdbr(species = "Homo sapiens", collection = "C5", subcollection = paste0("GO:", sub)),
    error = function(e) msigdbr(species = "Homo sapiens", category = "C5", subcategory = paste0("GO:", sub))
  )
  split(g$gene_symbol, g$gs_name)
}
pathways <- lapply(c(BP = "BP", CC = "CC", MF = "MF"), go_sets)

set.seed(42)
go_counts <- lapply(MAIN_CONTRASTS, function(cn) {
  ranks <- setNames(dep_df[[paste0("t_", cn)]], dep_df$gene)
  ranks <- ranks[!is.na(ranks) & !is.na(names(ranks)) & !duplicated(names(ranks))]
  ranks <- sort(ranks, decreasing = TRUE)
  lapply(names(ONT_COLORS), function(ont) {
    res <- tryCatch(
      fgsea(pathways[[ont]], ranks, minSize = 10, maxSize = 500),
      error = function(e) NULL
    )
    n <- if (is.null(res)) 0L else sum(res$padj < 0.05, na.rm = TRUE)
    tibble(contrast = cn, ontology = ont, n_terms = n)
  }) |> bind_rows()
}) |>
  bind_rows() |>
  mutate(
    contrast = factor(contrast, levels = MAIN_CONTRASTS),
    ontology = factor(ontology, levels = names(ONT_COLORS))
  )

write.csv(go_counts, file.path(DAT_DIR, "audit_panel_H_go_counts.csv"), row.names = FALSE)

totals <- go_counts |>
  group_by(contrast) |>
  summarise(tot = sum(n_terms), .groups = "drop")

PH_W <- 150
PH_H <- 120

# Per-contrast background bands tinted by the contrast code, grey borders between
band_df <- tibble(
  contrast = factor(MAIN_CONTRASTS, levels = MAIN_CONTRASTS),
  xmin = seq_along(MAIN_CONTRASTS) - 0.5,
  xmax = seq_along(MAIN_CONTRASTS) + 0.5,
  band = CONTRAST_COLORS[MAIN_CONTRASTS]
)

pH <- ggplot(go_counts, aes(contrast, n_terms, fill = ontology)) +
  geom_rect(
    data = band_df,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE, fill = scales::alpha(band_df$band, 0.14),
    color = "grey75", linewidth = 0.2
  ) +
  geom_col(width = 0.74, color = "black", linewidth = 0.25) +
  geom_text(
    data = totals, aes(contrast, tot, label = tot),
    inherit.aes = FALSE, vjust = -0.4, size = 2.8, fontface = "bold"
  ) +
  scale_fill_manual(values = ONT_COLORS) +
  scale_x_discrete(labels = CTR_SHORT) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
  labs(
    title = "Pathway Enrichment (GO GSEA)",
    subtitle = "Significant gene sets (BH < 0.05) by ontology",
    x = NULL, y = "Significant GO terms", fill = NULL, tag = "H"
  ) +
  FIG_THEME +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 30, hjust = 1))

save_panel(pH, file.path(RPT_DIR, "panel_h_pathways"), PH_W, PH_H)
cat("F02 Panel H done.\n")

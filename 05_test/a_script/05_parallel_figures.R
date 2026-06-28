# Parallel GSVA and ssGSEA figures: top divergent pathways per database, with the
# within-database FDR flagged and the nominal threshold marked. Same layout for
# both methods so they read side by side.
if (!exists("expr")) source(here::here("05_test", "a_script", "setup.R"))
res <- read_csv(file.path(DAT, "multidb_results.csv"), show_col_types = FALSE)

DB_ORDER <- c("Hallmark", "KEGG", "Reactome", "GO:BP", "GO Slim")

db_figure <- function(method_name, n_per_db = 6) {
  d <- res |>
    filter(method == method_name, !is.na(p_interaction)) |>
    mutate(database = factor(database, levels = DB_ORDER)) |>
    group_by(database) |>
    slice_min(p_interaction, n = n_per_db, with_ties = FALSE) |>
    ungroup() |>
    mutate(
      term = reorder_within(clean_pathway_name(pathway, 36), -log10(p_interaction), database),
      bias = if_else(slope > 0, "HR rises", "LR rises"),
      fdr_hit = fdr_within < 0.05
    )
  ggplot(d, aes(-log10(p_interaction), term, color = bias)) +
    geom_segment(aes(x = 0, xend = -log10(p_interaction), yend = term), linewidth = 0.5) +
    geom_point(aes(shape = fdr_hit), size = 2.2) +
    geom_vline(xintercept = -log10(0.05), linetype = "dotted", color = "grey55") +
    scale_color_manual(values = c("HR rises" = "#2166AC", "LR rises" = "#B2182B"), name = NULL) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 8), guide = "none") +
    scale_y_reordered() +
    facet_wrap(~database, scales = "free_y", ncol = 1) +
    labs(
      x = "-log10 p (group x timepoint)", y = NULL, title = method_name,
      subtitle = "star = within-database FDR < 0.05; dotted = nominal p = 0.05"
    ) +
    FIG_THEME +
    theme(axis.text.y = element_text(size = 6))
}

save_panel(db_figure("GSVA"), file.path(RPT, "fig5_gsva_by_db"), 180, 230)
save_panel(db_figure("ssGSEA"), file.path(RPT, "fig6_ssgsea_by_db"), 180, 230)

wide <- res |>
  select(pathway, database, method, p_interaction) |>
  pivot_wider(names_from = method, values_from = p_interaction)
hits <- c(
  "HALLMARK_HEME_METABOLISM", "HALLMARK_KRAS_SIGNALING_UP",
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE"
)
db_pal <- c(
  Hallmark = "#AA336A", KEGG = "#E65100", Reactome = "#1565C0",
  "GO:BP" = "#00796B", "GO Slim" = "#7B1FA2", Other = "grey60"
)
fig_agree <- ggplot(wide, aes(-log10(GSVA), -log10(ssGSEA), color = database)) +
  geom_abline(slope = 1, linetype = "dashed", color = "grey60") +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "grey70") +
  geom_vline(xintercept = -log10(0.05), linetype = "dotted", color = "grey70") +
  geom_point(size = 1.3, alpha = 0.55) +
  ggrepel::geom_text_repel(
    data = ~ filter(.x, pathway %in% hits),
    aes(label = clean_pathway_name(pathway, 24)), size = 2.4, color = "black",
    max.overlaps = 20, min.segment.length = 0
  ) +
  scale_color_manual(values = db_pal, name = NULL) +
  labs(
    x = "GSVA  -log10 p", y = "ssGSEA  -log10 p",
    title = "GSVA vs ssGSEA agreement on trajectory divergence",
    subtitle = sprintf("%d sets across 5 databases", nrow(wide))
  ) +
  FIG_THEME
save_panel(fig_agree, file.path(RPT, "fig4_method_agreement"), 185, 150)

cat("parallel figures written\n")

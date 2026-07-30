# Panel B: what the modules are. GO:BP over-representation of each module against the
# measured proteome as the universe. This is the one panel in F04 that carries a positive
# result, and it is a descriptive one: the modules are coherent biology, whatever they do or
# do not predict.
if (!exists("mods")) source(here::here("05_Figures", "shared", "WGCNA", "a_script", "setup.R"))
pacman::p_load(ggplot2, dplyr, forcats, stringr, shadowtext)

ora <- module_ora(mods, top_n = 3)

ora_plot <- ora |>
  filter(!is.na(term)) |>
  mutate(
    term = clean_pathway_name(term, max_chars = 40),
    module = fct_reorder(module, padj, .fun = min, .desc = TRUE),
    lp = -log10(padj)
  )

n_modules <- length(setdiff(unique(mods), "grey"))
n_with_term <- dplyr::n_distinct(ora_plot$module)

pB <- ggplot(ora_plot, aes(lp, reorder_within(term, lp, module), fill = module)) +
  geom_col(width = 0.75, color = "black", linewidth = 0.25) +
  shadowtext::geom_shadowtext(
    aes(label = term),
    x = 0.05, hjust = 0,
    size = 2.2, fontface = "bold", colour = "white", bg.colour = "grey20"
  ) +
  scale_fill_manual(values = setNames(levels(ora_plot$module), levels(ora_plot$module))) +
  scale_y_reordered() +
  facet_wrap(~module, scales = "free_y", ncol = 3) +
  labs(
    title = "What each module is (GO:BP over-representation)",
    subtitle = sprintf(
      "Universe = the measured proteome. Top 3 terms per module by BH q. %d of %d modules carry a significant term.",
      n_with_term, n_modules
    ),
    x = expression(-log[10] * " BH q"), y = NULL, tag = "B"
  ) +
  FIG_THEME +
  theme(
    legend.position = "none",
    axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.subtitle = element_text(size = 6.5, face = "italic", color = "grey40")
  )

save_panel(pB, file.path(RPT_DIR, "panels", "panel_b_pathways"), 230, 160)
F04_PANELS[["pathways"]] <- pB
F04_AUDIT[["module_ora"]] <- ora

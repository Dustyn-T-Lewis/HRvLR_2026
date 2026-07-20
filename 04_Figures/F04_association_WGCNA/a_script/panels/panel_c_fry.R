# Panel C: does any module move in any contrast?
#
# limma::fry, with the WGCNA modules as the index and the 03_DEP design, subject block and
# duplicateCorrelation carried through. This replaces a panel that ran fgsea with the modules
# as gene sets, ranked by moderated t from the matrix that defined them - circular, and it
# returned padj = 7e-33 in a study whose minimum q elsewhere is 0.36.
if (!exists("mods")) source(here::here("04_Figures", "F04_association_WGCNA", "a_script", "setup.R"))
pacman::p_load(ggplot2, dplyr, forcats, scales)

fry_res <- module_fry(
  mods, abund_dep, DESIGN, CONTRAST_MATRIX, BLOCK, CORRELATION
)

RESPONDER_CONTRASTS <- c(
  "Baseline_HRvLR", "Trained_HRvLR", "Acute_HRvLR",
  "Training_Interaction", "Acute_Interaction"
)

fry_plot <- fry_res |>
  mutate(
    family = if_else(contrast %in% RESPONDER_CONTRASTS, "Responder", "Within-group"),
    contrast = factor(contrast, levels = names(CTR_SHORT)),
    module = fct_reorder(module, fdr, .fun = min, .desc = TRUE),
    sig = fdr < 0.05
  )

pC <- ggplot(fry_plot, aes(contrast, module, fill = -log10(fdr))) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(
    aes(label = if_else(sig, "*", "")),
    size = 4, fontface = "bold", color = "grey15", vjust = 0.75
  ) +
  scale_fill_gradient(
    low = "grey93", high = "#2166AC", name = expression(-log[10] * " FDR"),
    limits = c(0, NA)
  ) +
  scale_x_discrete(labels = CTR_SHORT) +
  facet_grid(~family, scales = "free_x", space = "free_x") +
  labs(
    title = "Do the modules move? limma::fry, one test per module x contrast",
    subtitle = sprintf(
      "Rotation test on the 03_DEP design (subject-blocked, correlation = %.3f). * = FDR < 0.05. Minimum FDR across the responder contrasts: %.2f",
      CORRELATION,
      min(fry_res$fdr[fry_res$contrast %in% RESPONDER_CONTRASTS])
    ),
    x = NULL, y = NULL, tag = "C"
  ) +
  FIG_THEME +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.subtitle = element_text(size = 6.5, face = "italic", color = "grey40")
  )

save_panel(pC, file.path(RPT_DIR, "panels", "panel_c_fry"), 200, 130)
F04_PANELS[["fry"]] <- pC
F04_AUDIT[["module_fry"]] <- fry_res

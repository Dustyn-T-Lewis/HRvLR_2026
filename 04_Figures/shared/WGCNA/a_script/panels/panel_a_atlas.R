# Panel A: the module atlas. What the network found - module sizes, and how much of each
# module's eigengene variance is still explained by subject identity after the centred
# definition. The ICC column is the honest diagnostic: it is the number that exposed the
# earlier artifact, and it belongs on the face of the figure, not buried in a supplement.
#
# This panel owns two supplements: the soft-threshold sweep and the gene dendrogram, both of
# which document how the modules were arrived at.
if (!exists("mods")) source(here::here("04_Figures", "shared", "WGCNA", "a_script", "setup.R"))
pacman::p_load(ggplot2, dplyr, forcats, patchwork, scales, WGCNA)

# subject_variance keys modules by their ME* eigengene names; the colours are the bare form.
subj_var <- subject_variance(me_long) |>
  mutate(module = sub("^ME", "", module))

atlas <- tibble::enframe(mods, name = "gene", value = "module") |>
  filter(module != "grey") |>
  count(module, name = "n_proteins") |>
  left_join(subj_var, by = "module") |>
  mutate(module = fct_reorder(module, n_proteins))

p_size <- ggplot(atlas, aes(n_proteins, module, fill = module)) +
  geom_col(width = 0.75, color = "black", linewidth = 0.25) +
  geom_text(aes(label = n_proteins), hjust = -0.2, size = 2.4, fontface = "bold") +
  scale_fill_manual(values = setNames(levels(atlas$module), levels(atlas$module))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(title = "Module size", x = "Proteins", y = NULL, tag = "A") +
  FIG_THEME +
  theme(legend.position = "none", panel.grid.major.y = element_blank())

p_icc <- ggplot(atlas, aes(icc, module)) +
  geom_col(width = 0.75, fill = "grey75", color = "black", linewidth = 0.25) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  scale_x_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "Subject identity in the eigengene",
    subtitle = "ICC of ME by subject; dashed = half the variance",
    x = "ICC (subject)", y = NULL
  ) +
  FIG_THEME +
  theme(
    axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.subtitle = element_text(size = 6.5, face = "italic", color = "grey40")
  )

pA <- (p_size | p_icc) + plot_layout(widths = c(1, 1))

save_panel(pA, file.path(RPT_DIR, "panels", "panel_a_atlas"), 190, 110)
F04_PANELS[["atlas"]] <- pA
F04_AUDIT[["module_atlas"]] <- atlas

# Supplement: the soft-threshold sweep. Shows the fit clears RsquaredCut (so the FAQ's
# sample-size fallback table does not apply) and that mean connectivity at the chosen power
# is low, the second half of the Emory dual criterion.
sft <- wg$sft$fitIndices
sft$signed_r2 <- -sign(sft$slope) * sft$SFT.R.sq

p_r2 <- ggplot(sft, aes(Power, signed_r2)) +
  geom_hline(yintercept = 0.85, linetype = "dashed", color = "#B2182B", linewidth = 0.4) +
  geom_line(color = "grey60") +
  geom_point(aes(color = Power == wg$power), size = 2) +
  scale_color_manual(values = c(`TRUE` = "#2166AC", `FALSE` = "grey45"), guide = "none") +
  labs(
    title = "Scale-free fit", x = "Soft power", y = expression("signed " * R^2),
    subtitle = sprintf(
      "chosen beta = %d (first to clear RsquaredCut = 0.85); achieved R2 = %.3f",
      wg$power, wg$r2
    )
  ) +
  FIG_THEME +
  theme(plot.subtitle = element_text(size = 6.5, face = "italic", color = "grey40"))

p_k <- ggplot(sft, aes(Power, mean.k.)) +
  geom_line(color = "grey60") +
  geom_point(aes(color = Power == wg$power), size = 2) +
  scale_color_manual(values = c(`TRUE` = "#2166AC", `FALSE` = "grey45"), guide = "none") +
  scale_y_log10() +
  labs(
    title = "Mean connectivity", x = "Soft power", y = "mean k (log)",
    subtitle = sprintf("mean k at beta = %d is %.1f", wg$power, wg$mean_k)
  ) +
  FIG_THEME +
  theme(plot.subtitle = element_text(size = 6.5, face = "italic", color = "grey40"))

supp_sft <- (p_r2 | p_k)
save_panel(supp_sft, file.path(RPT_DIR, "supp", "panel_a_supp_soft_threshold"), 210, 100)
F04_AUDIT[["soft_threshold"]] <- sft[, c("Power", "signed_r2", "mean.k.", "median.k.", "max.k.")]

# Supplement: the gene dendrogram with the module colourbands beneath it.
png(
  file.path(RPT_DIR, "supp", "panel_a_supp_dendrogram.png"),
  width = 240, height = 130, units = "mm", res = 300, bg = "white"
)
WGCNA::plotDendroAndColors(
  wg$net$dendrograms[[1]],
  colors = wg$colors[wg$net$blockGenes[[1]]],
  groupLabels = "Module",
  dendroLabels = FALSE, hang = 0.03,
  addGuide = TRUE, guideHang = 0.05,
  main = "WGCNA gene dendrogram (within-subject-centred abundance)"
)
dev.off()

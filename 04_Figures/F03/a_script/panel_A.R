# F03 Panel A: WGCNA module map ----
setwd(rprojroot::find_rstudio_root_file())
if (!exists("wgcna_mem")) source("04_Figures/F03/a_script/HRvLR_F03_setup.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

# ── Module sizes (exclude grey) ----
mod_sizes <- wgcna_mem |>
  filter(group_id != "grey") |>
  count(module = group_id, name = "n_proteins") |>
  arrange(desc(n_proteins))

# Module order: descending protein count (top = most proteins on y-axis)
mod_order <- mod_sizes$module

# ── Phenotype associations (wgcna, tracks model) ----
trait_levels <- c(
  "comp_hypertrophy", "d_fcsa_I", "d_fcsa_II",
  "d_myovision_fcsa_I", "d_mcsa", "d_1rm_legpress", "d_1rm_ext"
)
trait_labels <- c(
  comp_hypertrophy    = "Composite hypertrophy",
  d_fcsa_I            = "D.fCSA I",
  d_fcsa_II           = "D.fCSA II",
  d_myovision_fcsa_I  = "D.MyoVision fCSA I",
  d_mcsa              = "D.mCSA",
  d_1rm_legpress      = "D.1RM leg-press",
  d_1rm_ext           = "D.1RM leg-ext"
)

assoc_df <- pheno_models |>
  filter(engine == "wgcna", model == "tracks") |>
  mutate(
    module = factor(module, levels = mod_order),
    trait_label = factor(trait_labels[trait], levels = unname(trait_labels))
  )

# Symmetric colour scale limit
beta_lim <- ceiling(max(abs(assoc_df$beta_std)) * 1000) / 1000
beta_lim <- max(beta_lim, 0.03)

# Audit CSV (long format)
audit_df <- assoc_df |>
  select(module, trait, trait_label, beta_std, p_bh)
write.csv(audit_df, file.path(DAT_DIR, "audit_panel_A_module_pheno.csv"), row.names = FALSE)

# ── Left panel: module-size bars ----
p_bars <- ggplot(mod_sizes, aes(x = n_proteins, y = factor(module, levels = mod_order), fill = module)) +
  geom_col(color = "black", linewidth = 0.3, width = 0.7) +
  geom_text(
    aes(label = n_proteins),
    hjust = -0.15, size = FIG_GEOM_TEXT, fontface = "bold"
  ) +
  scale_fill_identity() +
  scale_x_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(x = "Proteins (n)", y = NULL, title = NULL) +
  FIG_THEME +
  theme(
    axis.text.y = element_text(face = "bold"),
    panel.grid.major.y = element_blank()
  )

# ── Right panel: module × trait beta heatmap ----
# Bold x-axis label for the lead trait via label function
bold_if_lead <- function(x) {
  lapply(x, function(lbl) {
    if (lbl == "Composite hypertrophy") {
      bquote(bold(.(lbl)))
    } else {
      lbl
    }
  })
}

p_heat <- ggplot(assoc_df, aes(x = trait_label, y = module, fill = beta_std)) +
  geom_tile(color = "grey85", linewidth = 0.4) +
  geom_tile(
    data = filter(assoc_df, p_bh < 0.05),
    fill = NA, color = "black", linewidth = 1.0
  ) +
  geom_text(
    aes(label = sig_stars(p_bh)),
    size = 3, vjust = 0.75
  ) +
  scale_fill_gradient2(
    low = "#4393C3", mid = "white", high = "#D6604D",
    midpoint = 0,
    limits = c(-beta_lim, beta_lim),
    oob = scales::squish,
    name = "Std. b"
  ) +
  scale_x_discrete(position = "bottom", labels = bold_if_lead) +
  scale_y_discrete(limits = mod_order) +
  labs(x = NULL, y = NULL, title = NULL) +
  FIG_THEME +
  theme(
    axis.text.x = element_text(angle = 40, hjust = 1, size = FIG_AXIS_TEXT),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"
  )

# ── Compose ----
panel_a <- (p_bars | p_heat) +
  plot_layout(widths = c(1, 3)) +
  plot_annotation(
    title = "WGCNA proteome modules and their phenotype association",
    subtitle = sprintf(
      "n~16 subjects, hypothesis-generating; family-wise permutation p = %.3f (no module survives BH)",
      perm_null$p_perm[1]
    ),
    tag_levels = list(c("A", "")),
    theme = theme(
      plot.title    = element_text(face = "bold", size = FIG_TITLE_SIZE),
      plot.subtitle = element_text(face = "bold.italic", size = FIG_SUBTITLE_SIZE, color = "grey30"),
      plot.tag      = element_text(face = "bold", size = FIG_TAG_SIZE)
    )
  )

save_panel(panel_a, file.path(RPT_DIR, "panel_a_wgcna_map"), 210, 120)
cat("F03 Panel A done.\n")

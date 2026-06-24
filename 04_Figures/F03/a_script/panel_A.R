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

# Module order: descending protein count — largest (turquoise, 74) at top
mod_order <- mod_sizes$module

# ── Phenotype associations: all three wgcna models ----
trait_levels <- c(
  "comp_hypertrophy", "d_fcsa_I", "d_fcsa_II",
  "d_myovision_fcsa_I", "d_mcsa", "d_1rm_legpress", "d_1rm_ext"
)
trait_labels <- c(
  comp_hypertrophy    = "Composite hypertrophy",
  d_fcsa_I            = "fCSA I",
  d_fcsa_II           = "fCSA II",
  d_myovision_fcsa_I  = "MyoVision fCSA I",
  d_mcsa              = "mCSA",
  d_1rm_legpress      = "1RM leg-press",
  d_1rm_ext           = "1RM leg-ext"
)

model_labels <- c(
  tracks   = "Tracks (overall)",
  baseline = "Baseline -> response",
  acute    = "Acute shift"
)

assoc_df <- pheno_models |>
  filter(engine == "wgcna") |>
  mutate(
    module = factor(module, levels = rev(mod_order)),
    trait_label = factor(trait_labels[trait], levels = unname(trait_labels)),
    model = factor(model,
      levels = c("tracks", "baseline", "acute"),
      labels = unname(model_labels)
    )
  )

# Symmetric colour scale: driven by the widest-ranging model (acute, ~0.74)
wgcna_rows <- pheno_models[pheno_models$engine == "wgcna", ]
L <- ceiling(max(abs(wgcna_rows$beta_std)) * 10) / 10

# Significance subsets
thin_border <- filter(assoc_df, p < 0.05)
thick_border <- filter(assoc_df, p_bh < 0.05)

# Audit CSV (long format — all three models)
audit_df <- assoc_df |>
  select(engine = engine, module, trait, model, beta_std, p, p_bh) |>
  mutate(engine = "wgcna")

# Rebuild engine column from the source before mutate dropped it
audit_df <- pheno_models |>
  filter(engine == "wgcna") |>
  select(engine, module, trait, model, beta_std, p, p_bh)
write.csv(audit_df, file.path(DAT_DIR, "audit_panel_A_module_pheno.csv"), row.names = FALSE)

# ── Left panel: module-size bars ----
# Use rev(mod_order) so the largest module (turquoise) sits at the top
p_bars <- ggplot(
  mod_sizes,
  aes(
    x = n_proteins, y = factor(module, levels = rev(mod_order)),
    fill = module
  )
) +
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

# ── Right panel: module x trait heatmap, faceted by model ----
bold_if_lead <- function(x) {
  lapply(x, function(lbl) {
    if (lbl == "Composite hypertrophy") bquote(bold(.(lbl))) else lbl
  })
}

p_heat <- ggplot(assoc_df, aes(x = trait_label, y = module, fill = beta_std)) +
  geom_tile(color = "grey85", linewidth = 0.3) +
  geom_tile(
    data = thin_border,
    fill = NA, color = "grey20", linewidth = 0.4
  ) +
  geom_tile(
    data = thick_border,
    fill = NA, color = "black", linewidth = 1.1
  ) +
  geom_text(
    data = thick_border,
    aes(label = sig_stars(p_bh)),
    size = 3, vjust = 0.75
  ) +
  scale_fill_gradient2(
    low = "#4393C3", mid = "white", high = "#D6604D",
    midpoint = 0,
    limits = c(-L, L),
    oob = scales::squish,
    name = "Std. effect"
  ) +
  scale_x_discrete(position = "bottom", labels = bold_if_lead) +
  scale_y_discrete(limits = rev(mod_order)) +
  facet_wrap(~model, nrow = 1) +
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
  plot_layout(widths = c(1, 5)) +
  plot_annotation(
    title = "WGCNA proteome modules and their phenotype association",
    subtitle = sprintf(
      "n approx 16; traits = T1 to T2 change. Family-wise permutation p = %.3f; raw p<0.05 outlined, none survive BH.",
      perm_null$p_perm[1]
    ),
    tag_levels = list(c("A", "")),
    theme = theme(
      plot.title    = element_text(face = "bold", size = FIG_TITLE_SIZE),
      plot.subtitle = element_text(face = "bold.italic", size = FIG_SUBTITLE_SIZE, color = "grey30"),
      plot.tag      = element_text(face = "bold", size = FIG_TAG_SIZE)
    )
  )

save_panel(panel_a, file.path(RPT_DIR, "panel_a_wgcna_map"), 260, 120)
cat("F03 Panel A done.\n")

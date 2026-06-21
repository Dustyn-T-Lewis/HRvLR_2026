# F03 — Composite Assembly ----
# Sources all per-panel scripts, assembles with patchwork.
# Panels: A (DEP counts), B (Barcode), C (UpSet), D (fGSEA)
# Supp:   E (fGSEA stacked), S1.7 (pi-score distributions)

source("04_Figures/F03/a_script/panel_A.R")            # -> pA (+ sig_sets, dir_map for C)
source("04_Figures/F03/a_script/panel_B.R")            # -> pB
source("04_Figures/F03/a_script/panel_C_upset.R")      # -> pC_combined
source("04_Figures/F03/a_script/panel_D.R")            # -> pD
source("04_Figures/F03/a_script/supp/panel_E_stacked_fgsea.R")  # -> pE_combined (supp)
source("04_Figures/F03/a_script/supp/S1_7.R")          # -> s17 (supp)

message("All F03 panel scripts sourced — assembling composite...")

# Layout: left (A/D) | mid (C) | right (B)
left_col <- (pA / pD) +
  plot_layout(heights = c(0.30, 0.70))

fig <- (left_col | pC_combined | pB) +
  plot_layout(widths = c(0.28, 0.37, 0.35)) +
  plot_annotation(
    title = "Differential Expression & Pathway Enrichment",
    subtitle = "DEP counts, contrast overlap, barcode rank, and functional enrichment",
    theme = theme(
      plot.title    = element_text(face = "bold", size = 11, hjust = 0.5),
      plot.subtitle = element_text(size = 8, color = "grey30", hjust = 0.5)
    )
  )

ggsave(file.path(RPT_DIR, "F03_dep_pea_MAIN.pdf"), fig,
       width = 500, height = 280, units = "mm", device = PDF_DEVICE)
ggsave(file.path(RPT_DIR, "F03_dep_pea_MAIN.png"), fig,
       width = 500, height = 280, units = "mm", dpi = 300)
ggsave(file.path(RPT_DIR, "MAIN_F03_composite.pdf"), fig,
       width = 500, height = 280, units = "mm", device = PDF_DEVICE)
ggsave(file.path(RPT_DIR, "MAIN_F03_composite.png"), fig,
       width = 500, height = 280, units = "mm", dpi = 300)
cat("F03 composite saved\n")

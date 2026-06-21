# F02 — Composite Assembly ----
# Sources all per-panel scripts, assembles with patchwork.
# Panels: A (PCA), B (logFC), C (CV), D (CV comparison, supp), E (Imputed, supp)

source("04_Figures/F02/a_script/panel_A.R")     # -> pA_combined
source("04_Figures/F02/a_script/panel_B.R")     # -> pB
source("04_Figures/F02/a_script/panel_C.R")     # -> pC
source("04_Figures/F02/a_script/panel_D.R")     # -> pD (supp)
source("04_Figures/F02/a_script/panel_E.R")     # -> pE (supp)

message("All F02 panel scripts sourced — assembling composite...")

# Layout: top (A) | bottom-left (B) | bottom-right (C)
fig <- (pA_combined) /
  (pB | pC) +
  plot_layout(heights = c(0.40, 0.60)) +
  plot_annotation(
    title = "Proteomic Quality Control",
    subtitle = "PCA, effect size distributions, and inter-individual variability",
    theme = theme(
      plot.title    = element_text(face = "bold", size = 11, hjust = 0.5),
      plot.subtitle = element_text(size = 8, color = "grey30", hjust = 0.5)
    )
  )

ggsave(file.path(RPT_DIR, "F02_qc_MAIN.pdf"), fig,
       width = 300, height = 260, units = "mm", device = PDF_DEVICE)
ggsave(file.path(RPT_DIR, "F02_qc_MAIN.png"), fig,
       width = 300, height = 260, units = "mm", dpi = 300)
ggsave(file.path(RPT_DIR, "MAIN_F02_composite.pdf"), fig,
       width = 300, height = 260, units = "mm", device = PDF_DEVICE)
ggsave(file.path(RPT_DIR, "MAIN_F02_composite.png"), fig,
       width = 300, height = 260, units = "mm", dpi = 300)
cat("F02 composite saved\n")

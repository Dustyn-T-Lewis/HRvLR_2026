# F09 WGCNA — Figure-specific style
# Module-trait heatmap, eigengene triptych, hub networks
source("04_Figures/shared/style.R")
source("04_Figures/shared/pathway_utils.R")

# WGCNA data directories (produced by 04_Figures/F06/a_script/HRvLR_WGCNA_run.R)
WGCNA_DATA <- "04_Figures/F06/c_data/wgcna"
PANEL_DATA <- "04_Figures/F06/c_data"   # panel-ready intermediates

# Database colors for ORA bars
DB_COLORS <- c(

  "Hallmark"  = "#AA336A",
  "GO:BP"     = "#1565C0",
  "Reactome"  = "#00796B",
  "KEGG"      = "#E65100",
  "WikiPath"  = "#6A1B9A"
)

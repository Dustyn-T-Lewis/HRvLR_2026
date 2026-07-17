#!/usr/bin/env Rscript
# F03 build: the five responder rings as the main figure, the four within-group rings as
# the positive-control supplement, the per-contrast top-30 audits, and one workbook.
# The only script in this unit that writes.
pacman::p_load(here, patchwork, ggplot2, openxlsx, dplyr)

F03_PANELS <- list()
F03_REPORTS <- list()
F03_SUPP <- list()
F03_TABLES <- list()

a_script <- here("04_Figures", "F03_volcanoes", "a_script")
source(file.path(a_script, "setup.R"))
source(file.path(a_script, "rings.R"))
source(here("functions", "shared_utils.R"))

clear_dir(RPT_DIR)
supp_dir <- file.path(RPT_DIR, "supp")
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT_DIR, recursive = TRUE, showWarnings = FALSE)

source(file.path(a_script, "panels", "panel_a_responders.R"))
source(file.path(a_script, "panels", "panel_b_within_group.R"))
source(file.path(a_script, "composite.R"))

ggsave(file.path(RPT_DIR, "F03_volcanoes.png"), composite,
  width = 450, height = 320, units = "mm", dpi = 200, bg = "white"
)
ggsave(file.path(RPT_DIR, "F03_volcanoes.pdf"), composite,
  width = 450, height = 320, units = "mm", device = PDF_DEVICE, bg = "white"
)

save_panel(F03_PANELS$within_group, file.path(supp_dir, "within_group_responses"),
  width = 380, height = 320
)

for (grp in names(F03_SUPP)) {
  for (ct in names(F03_SUPP[[grp]])) {
    save_panel(F03_SUPP[[grp]][[ct]], file.path(supp_dir, paste0("top30_", ct)),
      width = 420, height = 250
    )
  }
}

fg_sig <- fg |>
  filter(!is.na(padj), padj < 0.05, is.finite(NES)) |>
  transmute(
    contrast, database,
    pathway = clean_pathway_name(pathway),
    NES, padj, size, leading_edge = leadingEdge
  ) |>
  arrange(contrast, padj)

ring_pathways <- bind_rows(F03_REPORTS) |>
  select(contrast, database, pathway, NES, padj, size, dedup_status,
    merged_into, overlap_jaccard, drawn,
    leading_edge = leadingEdge
  )

top30_updown <- bind_rows(F03_TABLES)

overview <- data.frame(
  sheet = c("fgsea_all", "fgsea_significant", "ring_pathways", "top30_updown"),
  description = c(
    "All fgsea enrichment across the 9 DEP contrasts, raw and un-deduplicated",
    "Significant subset (padj < 0.05), cleaned names, arranged by contrast and padj",
    "Ring-volcano pathways per contrast: EnrichmentMap dedup status, what merged, which arcs were drawn",
    "No-dedup top-30 up and top-30 down by NES for every contrast"
  ),
  stringsAsFactors = FALSE
)

wb <- createWorkbook()
for (s in c("overview", "fgsea_all", "fgsea_significant", "ring_pathways", "top30_updown")) {
  addWorksheet(wb, s)
}
writeData(wb, "overview", overview)
writeData(wb, "fgsea_all", fg)
writeData(wb, "fgsea_significant", fg_sig)
writeData(wb, "ring_pathways", ring_pathways)
writeData(wb, "top30_updown", top30_updown)
saveWorkbook(wb, file.path(DAT_DIR, "F03_volcanoes_source_data.xlsx"), overwrite = TRUE)

cat("F03 rebuilt: 5 responder rings, within-group supplement, top-30 audits, workbook\n")

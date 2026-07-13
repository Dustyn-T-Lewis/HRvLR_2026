#!/usr/bin/env Rscript
# F03 enrichVolcano ring-volcanoes. Build the fgsea cache (setup), render two
# composites -- the within-group responses beside the interactions, and the
# between-responder HRvLR differences -- plus one no-dedup top-30 figure per
# contrast, and append the display sheets to the source-data workbook.
pacman::p_load(here, patchwork, openxlsx, dplyr)

a_script <- here("04_Figures", "F03", "a_script")
source(file.path(a_script, "HRvLR_F03_setup.R"))
source(file.path(a_script, "ring_helpers.R"))
source(file.path(a_script, "panels", "panel_main.R"))
source(file.path(a_script, "panels", "panel_hrvlr.R"))
source(file.path(a_script, "panels", "supp", "panel_top30.R"))

supp_dir <- file.path(RPT_DIR, "supp")
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)

main <- panel_main(fg, dep, pw)
hrvlr <- panel_hrvlr(fg, dep, pw)

ggsave(file.path(RPT_DIR, "F03_composite.png"), main$plot,
  width = 540, height = 320, units = "mm", dpi = 200, bg = "white"
)
ggsave(file.path(RPT_DIR, "F03_composite.pdf"), main$plot,
  width = 540, height = 320, units = "mm", device = PDF_DEVICE, bg = "white"
)

ggsave(file.path(RPT_DIR, "F03_HRvLR.png"), hrvlr$plot,
  width = 450, height = 185, units = "mm", dpi = 200, bg = "white"
)
ggsave(file.path(RPT_DIR, "F03_HRvLR.pdf"), hrvlr$plot,
  width = 450, height = 185, units = "mm", device = PDF_DEVICE, bg = "white"
)

top30 <- panel_top30(fg)
for (ct in names(top30$plots)) {
  save_panel(top30$plots[[ct]], file.path(supp_dir, paste0("top30_", ct)),
    width = 420, height = 250
  )
}

ring_pathways <- bind_rows(main$reports, hrvlr$reports) |>
  select(contrast, database, pathway, NES, padj, size, dedup_status,
    merged_into, overlap_jaccard, drawn,
    leading_edge = leadingEdge
  )

overview <- data.frame(
  sheet = c("fgsea_all", "fgsea_significant", "ring_pathways", "top30_updown"),
  description = c(
    "All fgsea enrichment across the 9 DEP contrasts, raw and un-deduplicated; the cache F04/F05 read",
    "Significant subset (padj < 0.05), cleaned names, arranged by contrast and padj",
    "Ring-volcano pathways per contrast: EnrichmentMap dedup status, what merged, which arcs were drawn",
    "No-dedup top-30 up and top-30 down by FDR for every contrast"
  ),
  stringsAsFactors = FALSE
)

wb <- loadWorkbook(cache_xlsx)
for (s in c("overview", "ring_pathways", "top30_updown")) {
  if (s %in% names(wb)) removeWorksheet(wb, s)
}
addWorksheet(wb, "overview")
writeData(wb, "overview", overview)
addWorksheet(wb, "ring_pathways")
writeData(wb, "ring_pathways", ring_pathways)
addWorksheet(wb, "top30_updown")
writeData(wb, "top30_updown", top30$table)
worksheetOrder(wb) <- order(match(
  names(wb),
  c("overview", "fgsea_all", "fgsea_significant", "ring_pathways", "top30_updown")
))
saveWorkbook(wb, cache_xlsx, overwrite = TRUE)

message("F03 complete: 2 composites + per-contrast top-30 + workbook sheets.")

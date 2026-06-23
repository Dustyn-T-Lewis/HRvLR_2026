# F02 Panel G: contrast overlap of Pi < 0.05 DEPs (UpSet) ----
# Which proteins are shared across contrasts; mirrors YvO's UpSet overview.

setwd(rprojroot::find_rstudio_root_file())
if (!exists("meta")) source("04_Figures/F02/a_script/HRvLR_F02_setup.R")
suppressPackageStartupMessages(library(UpSetR))

dep_sets <- lapply(MAIN_CONTRASTS, function(c) {
  g <- dep_df$gene[which(dep_df[[paste0("sig_pi_", c)]] != 0)]
  unique(g[!is.na(g)])
})
names(dep_sets) <- unname(CTR_SHORT[MAIN_CONTRASTS])
dep_sets <- dep_sets[lengths(dep_sets) > 0]

file_stem <- file.path(RPT_DIR, "panel_g_upset")
draw_upset <- function() {
  print(upset(fromList(dep_sets),
    nsets = length(dep_sets), order.by = "freq", nintersects = 20,
    mainbar.y.label = "Shared DEPs", sets.x.label = "DEPs / contrast",
    main.bar.color = "grey25", sets.bar.color = "#4393C3",
    matrix.color = "#2166AC", text.scale = 1.1, point.size = 2.4
  ))
}
png(paste0(file_stem, ".png"), width = 160, height = 120, units = "mm", res = 300, bg = "white")
draw_upset()
dev.off()
pdf(paste0(file_stem, ".pdf"), width = 160 / 25.4, height = 120 / 25.4, bg = "white")
draw_upset()
dev.off()
cat("F02 Panel G done.\n")

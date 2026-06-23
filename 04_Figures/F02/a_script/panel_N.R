# F02 Panel N (pilot): 5-set overlap of the primary contrasts (Euler) ----
# Pi < 0.05 DEP sets for Baseline, Training HR/LR, Acute HR/LR; supplementary.

setwd(rprojroot::find_rstudio_root_file())
if (!exists("meta")) source("04_Figures/F02/a_script/HRvLR_F02_setup.R")
suppressPackageStartupMessages(library(eulerr))

primary <- c("Baseline_HRvLR", "Training_HR", "Training_LR", "Acute_HR", "Acute_LR")
sets <- lapply(primary, function(c) {
  g <- dep_df$gene[which(dep_df[[paste0("sig_pi_", c)]] != 0)]
  unique(g[!is.na(g)])
})
names(sets) <- unname(CTR_SHORT[primary])

fit <- euler(sets, shape = "ellipse")
p <- plot(fit,
  quantities = list(fontsize = 8),
  fills = list(fill = scales::alpha(unname(CONTRAST_COLORS[primary]), 0.40)),
  edges = list(col = unname(CONTRAST_COLORS[primary]), lwd = 1.5),
  labels = list(fontsize = 9, fontface = "bold"),
  main = list(label = "Primary contrasts (Pi < 0.05)", fontsize = 11, fontface = "bold")
)

stem <- file.path(RPT_DIR, "panel_n_venn_primary")
png(paste0(stem, ".png"), width = 110, height = 100, units = "mm", res = 300, bg = "white")
print(p)
dev.off()
pdf(paste0(stem, ".pdf"), width = 110 / 25.4, height = 100 / 25.4, bg = "white")
print(p)
dev.off()
cat("F02 Panel N done.\n")

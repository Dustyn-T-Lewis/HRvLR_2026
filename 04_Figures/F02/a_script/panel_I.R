# F02 Panel I: HR vs LR training-response DEP overlap (Venn) ----
# Shared vs responder-specific proteins changing with training (T2 - T1), Pi < 0.05.

setwd(rprojroot::find_rstudio_root_file())
if (!exists("meta")) source("04_Figures/F02/a_script/HRvLR_F02_setup.R")
suppressPackageStartupMessages(library(eulerr))

venn_panel <- function(hr_contrast, lr_contrast, title, file_stem,
                       w = 95, h = 95) {
  hr <- dep_df$gene[which(dep_df[[paste0("sig_pi_", hr_contrast)]] != 0)]
  lr <- dep_df$gene[which(dep_df[[paste0("sig_pi_", lr_contrast)]] != 0)]
  hr <- unique(hr[!is.na(hr)])
  lr <- unique(lr[!is.na(lr)])
  fit <- euler(list(HR = hr, LR = lr))
  p <- plot(fit,
    quantities = list(fontsize = 11),
    fills = list(fill = c(
      scales::alpha(GROUP_COLORS[["HR"]], 0.45),
      scales::alpha(GROUP_COLORS[["LR"]], 0.45)
    ), alpha = 1),
    edges = list(col = c(GROUP_COLORS[["HR"]], GROUP_COLORS[["LR"]]), lwd = 2),
    labels = list(fontsize = 12, fontface = "bold"),
    main = list(label = title, fontsize = 12, fontface = "bold")
  )
  for (dev_fun in list(
    function() png(paste0(file_stem, ".png"), width = w, height = h, units = "mm", res = 300, bg = "white"),
    function() pdf(paste0(file_stem, ".pdf"), width = w / 25.4, height = h / 25.4, bg = "white")
  )) {
    dev_fun()
    print(p)
    dev.off()
  }
  cat(sprintf(
    "  %s: HR=%d, LR=%d, shared=%d\n", basename(file_stem),
    length(hr), length(lr), length(intersect(hr, lr))
  ))
}

venn_panel(
  "Training_HR", "Training_LR",
  "Training Response (T2 - T1)",
  file.path(RPT_DIR, "panel_i_venn_training")
)
cat("F02 Panel I done.\n")

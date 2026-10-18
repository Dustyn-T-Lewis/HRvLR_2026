# F1 b: the proteome at baseline and its two responses, three PCAs in a row: T1 levels, each
# person's training change (T2 - T1) and acute change (T3 - T2). The T2 and T3 levels are in S1e.
library(here)
library(dplyr)
library(ggplot2)
library(patchwork)
source(here("04_Figures", "style", "figure_style.R"))
source(here("04_Figures", "style", "pca_grid.R"))

inputs <- pca_inputs
holds <- paste(
  "Panel b: PC1 and PC2 per person at T1 and for the training and acute changes, with each",
  "window's PERMANOVA for group, gap k and forced-split adjusted Rand index."
)
panel_data <- pca_data(c("T1", "Training", "Acute"))
time_test <- filter(read_sheet(inputs[["contrasts"]], "permanova"), term == "time")
stopifnot(nrow(time_test) == 1)

panel <- pca_row(
  panel_data,
  title = "Proteome at baseline and its response",
  subtitle = with(time_test, sprintf(
    paste0(
      "PCA of blood-corrected values (proteins scaled), per person; 95%% ellipse per group. ",
      "Per plot: PERMANOVA for HR/LR;\ngap k (1 = no natural clusters); ARI of a forced ",
      "two-cluster split against HR/LR (0 = chance). Time, all %d biopsies: R² %.2f, %s"
    ),
    n_samples, r2, label_p(p)
  ))
)

save_panel(panel, "F1", "b_pca", 155, 66)

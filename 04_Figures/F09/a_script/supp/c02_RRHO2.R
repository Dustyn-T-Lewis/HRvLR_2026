# F09 Supp C02: RRHO2 — Acute Concordance (Acute_HR vs Acute_LR)
# Stratified hypergeometric via RRHO2 package (Cahill et al. 2018), stepsize = 20
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F09/a_script/style.R")
source("04_Figures/shared/concordance_utils.R")
suppressPackageStartupMessages({
  library(tidyverse)
  library(RRHO2)
})

PG_W <- 180
RPT <- "04_Figures/F09/b_reports/supp"
RPT_PDF       <- file.path(RPT, "main", "pdf")
RPT_PNG       <- file.path(RPT, "main", "png")
RPT_SUPP_PDF  <- file.path(RPT, "supp", "pdf")
RPT_SUPP_PNG  <- file.path(RPT, "supp", "png")
DAT <- "04_Figures/F09/c_data"
dir.create(file.path(DAT, "panel_B"), recursive = TRUE, showWarnings = FALSE)
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF,      recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PNG,      recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_SUPP_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_SUPP_PNG, recursive = TRUE, showWarnings = FALSE)
pdf_device <- get_pdf_device()

JET_COLORS <- c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F",
                "yellow", "#FF7F00", "red", "#7F0000")

dep_df <- read_csv("03_DEP/a_non_imputed/c_data/03_combined_results.csv", show_col_types = FALSE)
rr_df <- dep_df %>%
  filter(!is.na(t_Acute_HR), !is.na(t_Acute_LR)) %>%
  distinct(gene, .keep_all = TRUE)

n_shared <- nrow(rr_df)

list1 <- data.frame(gene = rr_df$gene, score = rr_df$t_Acute_HR, stringsAsFactors = FALSE)
list2 <- data.frame(gene = rr_df$gene, score = rr_df$t_Acute_LR, stringsAsFactors = FALSE)

rrho_obj <- RRHO2_initialize(
  list1, list2,
  labels          = c("Acute_HR", "Acute_LR"),
  log10.ind       = TRUE,
  multipleTesting = "none",
  boundary        = 0.02,
  method          = "hyper",
  stepsize        = 20
)

hmat <- rrho_obj$hypermat
nr <- nrow(hmat); nc <- ncol(hmat)

grid_df <- expand.grid(i = seq_len(nr), j = seq_len(nc)) %>%
  mutate(neg_log10_p = as.vector(hmat))
max_val <- max(grid_df$neg_log10_p, na.rm = TRUE)

r_val   <- cor(rr_df$t_Acute_HR, rr_df$t_Acute_LR, method = "pearson")
rho_val <- cor(rr_df$t_Acute_HR, rr_df$t_Acute_LR, method = "spearman")
r_ci    <- fisher_z_ci(r_val, n_shared)

pB <- ggplot(grid_df, aes(x = i, y = j, fill = neg_log10_p)) +
  geom_raster() +
  scale_fill_gradientn(colors = JET_COLORS, limits = c(0, max_val),
                        na.value = "white",
                        name = expression(-log[10](p))) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_fixed() +
  labs(title = "RRHO2: Acute Concordance (HR vs LR)",
       subtitle = sprintf("Stratified hypergeometric (stepsize=20) | N = %d | r = %.3f [%.3f, %.3f] | rho = %.3f",
                            n_shared, r_val, r_ci[1], r_ci[2], rho_val),
       x = "Acute HR rank (Up -> Down)",
       y = "Acute LR rank (Up -> Down)") +
  FIG_THEME +
  theme(legend.position = "right", legend.key.width = unit(3, "mm"),
        legend.key.height = unit(15, "mm"),
        axis.text = element_blank(), axis.ticks = element_blank())

ggsave(file.path(RPT_SUPP_PDF, "c02_RRHO2_SUPP.pdf"), pB,
       width = PG_W, height = PG_W, units = "mm", device = pdf_device)
ggsave(file.path(RPT_SUPP_PNG, "c02_RRHO2_SUPP.png"), pB,
       width = PG_W, height = PG_W, units = "mm", dpi = 300)

grid_df %>%
  transmute(grid_x = i, grid_y = j, neg_log10_p = round(neg_log10_p, 3)) %>%
  write_csv(file.path(DAT, "panel_B", "rrho2_grid.csv"))

message("F09 Supp C02 RRHO2 done")

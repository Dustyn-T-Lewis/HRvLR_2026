# F06 Panel B (HR): RRHO2 — Chronic vs Acute Ranked Overlap
# Training_HR (x) vs Acute_HR (y)
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F06/a_script/style.R")

suppressPackageStartupMessages({ library(tidyverse) })

GRP <- "HR"; CTR_X <- "Training_HR"; CTR_Y <- "Acute_HR"
PG_W <- 180
RPT <- "04_Figures/F06/b_reports"
DAT <- "04_Figures/F06/c_data"
dir.create(file.path(DAT, paste0("panel_B_", GRP)), recursive = TRUE, showWarnings = FALSE)
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
pdf_device <- get_pdf_device()

dep_df <- read_csv("03_DEP/a_non_imputed/c_data/03_combined_results.csv", show_col_types = FALSE)

rank_df <- dep_df %>%
  filter(!is.na(.data[[paste0("t_", CTR_X)]]), !is.na(.data[[paste0("t_", CTR_Y)]])) %>%
  mutate(rank_X = rank(-.data[[paste0("t_", CTR_X)]], ties.method = "average"),
         rank_Y = rank(-.data[[paste0("t_", CTR_Y)]], ties.method = "average"))

rrho <- compute_rrho2_grid(rank_df$rank_X, rank_df$rank_Y)
mat <- rrho$matrix; brk <- rrho$breaks; n <- rrho$n

grid_df <- expand.grid(i = seq_len(nrow(mat)), j = seq_len(ncol(mat))) %>%
  mutate(x = brk[i], y = brk[j], signed_score = mat[cbind(i, j)])
max_abs <- max(abs(grid_df$signed_score))

r_val  <- cor(rank_df[[paste0("t_", CTR_X)]], rank_df[[paste0("t_", CTR_Y)]], method = "pearson")
rho_val <- cor(rank_df[[paste0("t_", CTR_X)]], rank_df[[paste0("t_", CTR_Y)]], method = "spearman")
r_ci <- fisher_z_ci(r_val, nrow(rank_df))

pB <- ggplot(grid_df, aes(x = x, y = y, fill = signed_score)) +
  geom_tile() +
  scale_fill_gradient2(low = "#4393C3", mid = "white", high = "#D6604D",
                        midpoint = 0, limits = c(-max_abs, max_abs),
                        name = expression(-log[10](p) %*% sign)) +
  geom_hline(yintercept = n / 2, linetype = "dashed", color = "grey40", linewidth = 0.3) +
  geom_vline(xintercept = n / 2, linetype = "dashed", color = "grey40", linewidth = 0.3) +
  annotate("text", x = n * 0.25, y = n * 0.75, label = "Concordant\nUp",
           fontface = "bold", size = 4, color = "#D6604D") +
  annotate("text", x = n * 0.75, y = n * 0.25, label = "Concordant\nDown",
           fontface = "bold", size = 4, color = "#D6604D") +
  annotate("text", x = n * 0.75, y = n * 0.75, label = "Discordant",
           fontface = "bold", size = 4, color = "#4393C3") +
  annotate("text", x = n * 0.25, y = n * 0.25, label = "Discordant",
           fontface = "bold", size = 4, color = "#4393C3") +
  coord_fixed() +
  labs(title = sprintf("RRHO2: Chronic vs Acute (%s)", GRP),
       subtitle = sprintf("N = %d | r = %.3f [%.3f, %.3f] | rho = %.3f",
                            n, r_val, r_ci[1], r_ci[2], rho_val),
       x = paste0("Rank by t(", CTR_FACET[CTR_X], ") -> most up"),
       y = paste0("Rank by t(", CTR_FACET[CTR_Y], ") -> most up")) +
  FIG_THEME +
  theme(legend.position = "right", legend.key.width = unit(3, "mm"),
        legend.key.height = unit(15, "mm"), axis.text = element_text(size = 8))

ggsave(file.path(RPT, paste0("panel_B_RRHO2_", GRP, "_MAIN.pdf")), pB,
       width = PG_W, height = PG_W, units = "mm", device = pdf_device)
ggsave(file.path(RPT, paste0("panel_B_RRHO2_", GRP, "_MAIN.png")), pB,
       width = PG_W, height = PG_W, units = "mm", dpi = 300)

grid_df %>%
  transmute(rank_X = x, rank_Y = y, signed_log10p = round(signed_score, 3)) %>%
  write_csv(file.path(DAT, paste0("panel_B_", GRP), "rrho2_grid.csv"))

message(sprintf("F06 Panel B RRHO2 (%s) done", GRP))

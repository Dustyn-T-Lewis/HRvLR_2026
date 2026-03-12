# F0 Panel A: Accumulated Volume Load ----
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggsignif)
})

PW <- 85
PH <- 80

RPT <- "04_Figures/F0/b_reports"
DAT <- "04_Figures/F0/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

# Data ----
meta <- read_csv("00_input/HRvLR_meta.csv", show_col_types = FALSE)

vl_df <- meta %>%
  filter(Timepoint == "T3", !is.na(ACCUM_VL)) %>%
  select(Subject_ID, Group, ACCUM_VL) %>%
  mutate(Group = factor(Group, levels = c("HR", "LR")))

stats_A <- t.test(ACCUM_VL ~ Group, data = vl_df)

# Plot ----
bar_colors <- c(HR = unname(GROUP_FILL["HR_T2"]),
                LR = unname(GROUP_FILL["LR_T2"]))

pA <- ggplot(vl_df, aes(x = Group, y = ACCUM_VL, fill = Group)) +
  geom_bar(stat = "summary", fun = mean, width = 0.45,
           color = "black", linewidth = 0.3) +
  geom_errorbar(stat = "summary", fun.data = mean_se,
                width = 0.2, linewidth = 0.4) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.5,
              shape = 16, color = "grey30") +
  geom_signif(
    comparisons = list(c("HR", "LR")),
    annotations = fmt_p(stats_A$p.value),
    textsize = KEY_TEXT, tip_length = 0.02,
    y_position = bracket_pos(vl_df$ACCUM_VL)
  ) +
  scale_fill_manual(values = bar_colors) +
  scale_x_discrete(labels = c(
    HR = sprintf("HR (n = %d)", sum(vl_df$Group == "HR")),
    LR = sprintf("LR (n = %d)", sum(vl_df$Group == "LR")))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)),
                     labels = scales::label_comma()) +
  labs(title = "Accumulated Volume Load", y = "Volume load (kg)", x = NULL,
       tag = "A") +
  coord_cartesian(clip = "off") +
  FIG_THEME + theme(legend.position = "none",
                    plot.margin = margin(5, 5, 5, 5))

# Audit ----
audit_A <- vl_df %>%
  group_by(Group) %>%
  summarise(n = n(), mean = mean(ACCUM_VL), sd = sd(ACCUM_VL),
            sem = sd(ACCUM_VL) / sqrt(n()), .groups = "drop") %>%
  mutate(t_test_p = stats_A$p.value)
write_csv(audit_A, file.path(DAT, "audit_panel_A_volume_load.csv"))

# Save ----
save_panel(pA, file.path(RPT, "panel_A_volume_load"), PW, PH)

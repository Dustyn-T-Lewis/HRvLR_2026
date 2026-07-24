# Assemble the F04 association composite from the leaf workbooks: the variance
# partition, the global ordination test, and null-summary grids showing that no
# association survives BH in any contrast x method x feature space. In-sample
# effects exist and are shown continuously; the grids report how many clear BH
# (expected zero at n=16) and the smallest q reached.

pacman::p_load(here, dplyr, purrr, tidyr, ggplot2, openxlsx, patchwork)
source(here("functions", "shared_style.R"))
source(here("functions", "assoc_panels.R"))

root <- here("04_Figures", "F04_association")

read_sheet <- function(path, sheet) {
  if (!file.exists(path) || !sheet %in% getSheetNames(path)) {
    return(NULL)
  }
  read.xlsx(path, sheet = sheet)
}

# variance partition
vp <- read_sheet(
  file.path(root, "global", "mixed_model", "c_data", "source_data.xlsx"),
  "variance_partition"
)
p_vp <- build_varpart_panel(vp)

# PERMANOVA line for the subtitle
perm <- read_sheet(
  file.path(root, "global", "ordination", "c_data", "source_data.xlsx"),
  "permanova"
)
perm_group <- perm |> filter(term == "Group")
ord_line <- sprintf(
  "PERMANOVA group R^2 = %.3f, p = %.2f (strata = subject)",
  perm_group$R2, perm_group$`Pr(>F)`
)

# group DE null grid
group_null <- map_dfr(c("T1", "T2", "T3", "training", "acute"), function(ct) {
  map_dfr(c("limma", "mixed_slice"), function(mt) {
    wb <- file.path(root, "group_HRvLR", ct, mt, "c_data", "source_data.xlsx")
    map_dfr(c("proteins", "pathways", "modules"), function(sp) {
      d <- read_sheet(wb, sp)
      if (is.null(d)) {
        return(NULL)
      }
      tibble(
        contrast = ct, method = mt, space = sp,
        n_hit = sum(d$bh < 0.05, na.rm = TRUE),
        min_bh = min(d$bh, na.rm = TRUE)
      )
    })
  })
}) |>
  mutate(
    contrast = factor(
      contrast,
      levels = c("T1", "T2", "T3", "training", "acute")
    ),
    space = factor(space, levels = c("proteins", "pathways", "modules")),
    row = interaction(method, space, sep = " / "),
    lab = sprintf("%d\n(q=%.2f)", n_hit, min_bh)
  )

p_group <- ggplot(group_null, aes(contrast, row, fill = n_hit)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = lab), size = 2.2, lineheight = 0.8) +
  scale_fill_gradient(low = "grey92", high = "#B2182B", name = "BH hits") +
  labs(
    x = NULL, y = NULL, title = "Group HR vs LR — BH survivors per contrast",
    subtitle = "Count clearing BH q<.05 and the smallest q reached."
  ) +
  FIG_THEME +
  theme(
    plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE),
    axis.text.y = element_text(size = 7)
  )

# phenotype null grid (lead trait: composite hypertrophy)
pheno_null <- map_dfr(c("baseline", "training", "acute"), function(ph) {
  wb <- file.path(root, "phenotype", ph, "limma", "c_data", "source_data.xlsx")
  map_dfr(c("proteins", "pathways", "modules"), function(sp) {
    d <- read_sheet(wb, paste0("assoc_", sp))
    if (is.null(d)) {
      return(NULL)
    }
    d <- d |> filter(trait == "comp_hypertrophy")
    tibble(
      phase = ph, space = sp,
      n_hit = sum(d$bh < 0.05, na.rm = TRUE),
      min_bh = min(d$bh, na.rm = TRUE)
    )
  })
}) |>
  mutate(
    phase = factor(phase, levels = c("baseline", "training", "acute")),
    space = factor(space, levels = c("proteins", "pathways", "modules")),
    lab = sprintf("%d\n(q=%.2f)", n_hit, min_bh)
  )

p_pheno <- ggplot(pheno_null, aes(phase, space, fill = n_hit)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = lab), size = 2.3, lineheight = 0.8) +
  scale_fill_gradient(low = "grey92", high = "#238B45", name = "BH hits") +
  labs(
    x = NULL, y = NULL,
    title = "Phenotype (vs composite hypertrophy) — BH survivors",
    subtitle = "Acute reflects damage-repair, not hypertrophy (Damas 2016)."
  ) +
  FIG_THEME +
  theme(plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE))

composite <- (p_vp | p_group) / p_pheno +
  plot_layout(heights = c(1, 0.7)) +
  plot_annotation(
    title = "F04 — proteome association with training response (in-sample)",
    subtitle = paste0(
      "Descriptive; F05 adjudicates generalization. ", ord_line, "."
    ),
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(face = "bold", size = FIG_TITLE_SIZE + 1),
      plot.subtitle = element_text(face = "italic", colour = "grey30")
    )
  )

dir.create(file.path(root, "b_reports"), showWarnings = FALSE)
dir.create(file.path(root, "c_data"), showWarnings = FALSE)
save_panel(composite, file.path(root, "b_reports", "F04_composite"),
  width = 250, height = 240
)

wb <- createWorkbook()
addWorksheet(wb, "group_null")
writeData(wb, "group_null", group_null |> select(-row, -lab))
addWorksheet(wb, "phenotype_null")
writeData(wb, "phenotype_null", pheno_null |> select(-lab))
addWorksheet(wb, "variance_partition")
writeData(wb, "variance_partition", vp)
saveWorkbook(wb, file.path(root, "c_data", "F04_metrics.xlsx"),
  overwrite = TRUE
)

message("F04 group BH survivors total: ", sum(group_null$n_hit))
message("F04 phenotype BH survivors total: ", sum(pheno_null$n_hit))

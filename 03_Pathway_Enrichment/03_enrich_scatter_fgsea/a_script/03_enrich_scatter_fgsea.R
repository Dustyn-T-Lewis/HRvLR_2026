# Put HR and LR on one pair of axes, once for training and once for the acute bout. fgsea scores
# each set once per contrast, so one arm's NES against the other's compares the arms set by set.
# No new test; it reshapes 01_run_fgsea_and_fry output.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(ggplot2)
  library(patchwork)
})

stage <- here("03_Pathway_Enrichment", "03_enrich_scatter_fgsea")
figure_dir <- file.path(stage, "b_reports")
out <- file.path(stage, "c_data")
walk(c(figure_dir, out), dir.create, recursive = TRUE, showWarnings = FALSE)
# Clear last run's figures so the bundle holds only this run's pages.
unlink(list.files(figure_dir, "[.](png|pdf)$", full.names = TRUE))

inputs <- c(set_tests = "03_Pathway_Enrichment/01_run_fgsea_and_fry/c_data/set_tests.rds")
paths <- map_chr(inputs, here)
if (!all(file.exists(paths))) {
  stop("Run 01_run_fgsea_and_fry first. Missing: ", inputs[["set_tests"]])
}
fg <- readRDS(paths[["set_tests"]])
manifest <- tibble(
  input = names(inputs), path = unname(inputs), md5 = unname(tools::md5sum(paths))
)

run_pair <- function(x_contrast, y_contrast, tag, prefix) {
  paired <- fg$set_tests |>
    filter(method == "fgsea", contrast %in% c(x_contrast, y_contrast)) |>
    mutate(contrast = if_else(contrast == x_contrast, "x", "y")) |>
    pivot_wider(
      id_cols = c(set_id, database, pathway),
      names_from = contrast, values_from = c(n, NES, padj, main)
    ) |>
    mutate(
      # A protein untested in one contrast leaves that ranking, so a set's measured size can
      # differ by one or two between arms; the point is sized by the larger.
      n = pmax(n_x, n_y),
      significance = case_when(
        padj_x < 0.05 & padj_y < 0.05 ~ "Both",
        padj_x < 0.05 ~ x_contrast,
        padj_y < 0.05 ~ y_contrast,
        .default = "NS"
      ) |> factor(c("Both", x_contrast, y_contrast, "NS")),
      # Opposite signs in the two arms, whichever arm, if either, is significant.
      discordant = sign(NES_x) != sign(NES_y),
      survivor = main_x | main_y,
      label = enrichVolcano::ev_clean_label(pathway)
    )
  # A set near the size floor can drop out of one arm entirely. Only sets scored in both are
  # compared.
  paired <- filter(paired, !is.na(NES_x), !is.na(NES_y))
  stopifnot(nrow(paired) > 0)

  concordance <- function(data) {
    # A quadrant panel can hold one or two sets, too few for cor.test.
    test <- if (nrow(data) >= 3) {
      cor.test(data$NES_x, data$NES_y, method = "spearman", exact = FALSE)
    } else {
      list(estimate = NA_real_, p.value = NA_real_)
    }
    hits <- filter(data, significance != "NS")
    tibble(
      sets = nrow(data), significant = nrow(hits), discordant = sum(hits$discordant),
      rho = round(unname(test$estimate), 3), p = test$p.value,
      same_sign = round(mean(!hits$discordant), 3)
    )
  }

  set_colours <- set_names(
    c("#6A3D9A", "#D7301F", "#2B6CB0"), c("Both", x_contrast, y_contrast)
  )

  # One panel builder for all six panels. `labelled` is the subset that gets names, so a dense
  # cloud and a zoomed handful of sets differ only in what is passed in.
  nes_panel <- function(data, title, labelled = data[0, ], pad = 0.12) {
    span <- range(c(data$NES_x, data$NES_y))
    limit <- span + c(-1, 1) * diff(span) * pad
    shown <- filter(data, significance != "NS")
    stats <- concordance(data)
    # A panel reports rho only with 30 or more sets.
    caption <- sprintf(
      "%d set%s | %d significant%s | %d discordant", stats$sets, if (stats$sets == 1) "" else "s",
      stats$significant,
      if (stats$sets >= 30) sprintf(" | rho %.2f", stats$rho) else "", stats$discordant
    )
    ggplot(data, aes(NES_x, NES_y)) +
      geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
      geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
      geom_abline(slope = 1, linetype = "dashed", colour = "grey45", linewidth = 0.4) +
      geom_point(
        data = filter(data, significance == "NS"),
        colour = "grey82", size = 0.45, alpha = 0.3
      ) +
      geom_point(aes(colour = significance, size = n), data = shown, alpha = 0.85) +
      ggrepel::geom_text_repel(
        data = labelled, aes(label = label), size = 2.3, colour = "grey15",
        segment.colour = "grey60", segment.size = 0.25, min.segment.length = 0,
        max.overlaps = Inf, seed = 1, force = 6, lineheight = 0.85
      ) +
      scale_colour_manual(values = set_colours, name = "significant in", drop = FALSE) +
      scale_size_continuous(
        range = c(1, 4.5), name = "genes",
        limits = range(paired$n), breaks = c(50, 150, 300)
      ) +
      coord_fixed(xlim = limit, ylim = limit) +
      labs(
        x = paste("NES,", x_contrast), y = paste("NES,", y_contrast),
        title = title, subtitle = caption
      ) +
      theme_minimal(base_size = 9) +
      theme(
        plot.title = element_text(face = "bold", size = 10),
        plot.subtitle = element_text(size = 7.5, colour = "grey30")
      )
  }

  save_composite <- function(figure, name, width, height) {
    walk(c("png", "pdf"), \(extension) {
      ggsave(file.path(figure_dir, paste0(name, ".", extension)), figure,
        width = width, height = height, dpi = 300, bg = "white"
      )
    })
    message("wrote ", name)
  }

  page_labels <- function(title, subtitle, caption) {
    plot_annotation(
      title = title, subtitle = subtitle, caption = caption, tag_levels = "A",
      theme = theme(
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 8.5, colour = "grey30"),
        plot.caption = element_text(size = 7.5, colour = "grey45", hjust = 0)
      )
    )
  }

  # Composite one: every collection, then the sets that survived collapse, then the discordant
  # handful on their own axes.
  discordant_sets <- filter(paired, discordant, significance != "NS")
  survivors <- filter(paired, survivor)
  composite_all <- wrap_plots(
    nes_panel(paired, "All collections"),
    nes_panel(survivors, "Collapse survivors"),
    # Too few points for a size key, and patchwork will not merge guide sets that differ.
    nes_panel(discordant_sets, "Discordant", labelled = discordant_sets, pad = 0.35) +
      guides(colour = "none", size = "none"),
    nrow = 1
  ) +
    page_labels(
      paste0("NES concordance, all collections: ", x_contrast, " against ", y_contrast),
      sprintf("fgsea NES per contrast, %d sets, BH within contrast", nrow(paired)),
      paste(
        "One point per gene set, scored in both", tag, "contrasts. Dashed line is identity;",
        "grey points reach neither threshold. Discordant sets are significant in at least one",
        "arm and sit on opposite sides of zero. Panel C rescales them.",
        "Table: c_data/nes_scatter.csv."
      )
    ) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  save_composite(composite_all, paste0(prefix, "_nes_concordance_all_", tag), 13, 6)

  # Composite two: the two collections whose members do not nest, then each concordant quadrant
  # scaled to its own points so every set can be named.
  curated <- filter(paired, database %in% c("Hallmark", "GO_Slim"))
  quadrant <- function(direction) {
    rows <- filter(curated, significance != "NS", !discordant, (NES_x > 0) == direction)
    # Panel A carries the colour key; a quadrant missing a category would draw a second one.
    nes_panel(rows, if (direction) "Up in both" else "Down in both", labelled = rows, pad = 0.22) +
      guides(colour = "none")
  }
  composite_curated <- (
    nes_panel(
      curated, "Hallmark and GO Slim",
      labelled = slice_min(
        filter(curated, significance != "NS"), padj_x + padj_y,
        n = 8, with_ties = FALSE
      )
    ) | (quadrant(TRUE) / quadrant(FALSE))
  ) +
    page_labels(
      paste0("NES concordance, Hallmark and GO Slim: ", x_contrast, " against ", y_contrast),
      sprintf("fgsea NES per contrast, %d non-nesting sets, BH within contrast", nrow(curated)),
      paste(
        "Panel A is every Hallmark and GO Slim set; B and C rescale the significant sets with the",
        "same sign in both arms, so each can be named. Point size is gene count, colour the",
        "contrast a set reached FDR 0.05 in. Table: c_data/nes_scatter.csv."
      )
    ) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  save_composite(composite_curated, paste0(prefix, "_nes_concordance_curated_", tag), 12, 7)

  summary_table <- bind_rows(
    mutate(concordance(paired), population = "all collections"),
    mutate(concordance(survivors), population = "collapse survivors"),
    mutate(concordance(curated), population = "Hallmark and GO Slim")
  ) |>
    relocate(population)
  export <- paired |>
    transmute(
      pair = tag, x_contrast, y_contrast, set_id, database, pathway, label,
      genes = n,
      nes_x = round(NES_x, 3), nes_y = round(NES_y, 3),
      padj_x = signif(padj_x, 4), padj_y = signif(padj_y, 4),
      significance = as.character(significance), discordant, survivor
    ) |>
    arrange(significance, desc(abs(nes_x) + abs(nes_y)))
  list(summary = mutate(summary_table, pair = tag, .before = 1), export = export)
}

pairs <- list(
  run_pair("Training_HR", "Training_LR", "training", "01"),
  run_pair("Acute_HR", "Acute_LR", "acute", "02")
)
summary_table <- list_rbind(map(pairs, "summary"))
export <- list_rbind(map(pairs, "export"))
print(as.data.frame(summary_table))

packages <- c("here", "dplyr", "tidyr", "ggplot2", "ggrepel", "patchwork", "enrichVolcano")
versions <- tibble(
  package = packages, version = map_chr(packages, \(p) as.character(packageVersion(p)))
)
readr::write_csv(export, file.path(out, "nes_scatter.csv"))
writexl::write_xlsx(
  list(
    concordance = summary_table,
    discordant = filter(export, discordant, significance != "NS"),
    nes_scatter = export, input_manifest = manifest, package_versions = versions
  ),
  file.path(out, "03_enrich_scatter_fgsea.xlsx")
)
combined <- file.path(figure_dir, "03_enrich_scatter_fgsea_figures.pdf")
pages <- setdiff(list.files(figure_dir, "[.]pdf$", full.names = TRUE), combined)
invisible(qpdf::pdf_combine(sort(pages), combined))
message(
  "wrote 03_enrich_scatter_fgsea.xlsx and ", length(pages),
  " figures bundled into ", qpdf::pdf_length(combined), " pages"
)

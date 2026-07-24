# Leaf drivers for the F04 association suite. The global mixed model is fit once
# and cached; every group_HRvLR mixed-slice leaf reads that cache, so the five
# group contrasts are slices of one fit rather than five independent models
# (methods doc Part 8). limma leaves add the leaner per-contrast snapshot and
# reconcile against 03_DEP.

pacman::p_load(here, dplyr, tibble, openxlsx, ggplot2, patchwork)
source(here("functions", "shared_style.R"))
source(here("functions", "assoc_features.R"))
source(here("functions", "assoc_panels.R"))

ASSOC_SPACES <- c("proteins", "pathways", "modules")
PHENO_TRAITS <- c(
  "comp_hypertrophy", "d_fcsa_I", "d_fcsa_II", "d_mcsa",
  "d_1rm_legpress", "d_1rm_ext"
)

global_cache_path <- function() {
  here(
    "04_Figures", "F04_association", "global", "mixed_model", "c_data",
    "global_contrasts.rds"
  )
}

leaf_dirs <- function(out_dir) {
  dat <- file.path(out_dir, "c_data")
  rpt <- file.path(out_dir, "b_reports")
  dir.create(dat, recursive = TRUE, showWarnings = FALSE)
  dir.create(rpt, recursive = TRUE, showWarnings = FALSE)
  list(dat = dat, rpt = rpt)
}

write_leaf_workbook <- function(path, sheets) {
  wb <- createWorkbook()
  for (nm in names(sheets)) {
    addWorksheet(wb, nm)
    writeData(wb, nm, sheets[[nm]])
  }
  saveWorkbook(wb, path, overwrite = TRUE)
}

# Global mixed model: fit `~ group*timepoint + (1|subject)` once per feature in
# each space, cache the six contrasts for the group leaves, and pair the protein
# variance partition with the group-main and training-interaction volcanoes.
run_global_mixed_leaf <- function(bundle, out_dir, cores = 2L) {
  dirs <- leaf_dirs(out_dir)
  contrasts_by_space <- lapply(bundle$feature_sets, function(mat) {
    associate_global_hlm(mat, bundle$meta, cores = cores)
  })
  saveRDS(contrasts_by_space, global_cache_path())

  vp <- varpart_global(bundle$imputed_proteins, bundle$meta)
  prot <- contrasts_by_space$proteins

  slice_space <- function(space, contrast) {
    contrasts_by_space[[space]] |>
      filter(.data$contrast == !!contrast) |>
      transmute(feature, effect = estimate, t, p, bh)
  }

  p_vp <- build_varpart_panel(vp)
  p_gm <- assoc_volcano(slice_space("proteins", "group_main"), bundle$gene_map,
    title = "Group main effect (proteins)"
  )
  p_int <- assoc_volcano(slice_space("proteins", "training"), bundle$gene_map,
    title = "Training interaction (proteins)",
    effect_lab = "interaction (HR - LR)"
  )
  p_pw <- assoc_volcano(slice_space("pathways", "group_main"),
    title = "Group main effect (pathways)", effect_lab = "effect (HR - LR)"
  )

  panel <- (p_vp | p_gm) / (p_int | p_pw) +
    plot_annotation(
      title = "Global association — one mixed model, all samples",
      subtitle = paste(
        "Descriptive, in-sample. Fixed effects trusted; variance",
        "components descriptive (n=16). F05 adjudicates."
      ),
      theme = theme(
        plot.title = element_text(face = "bold", size = FIG_TITLE_SIZE),
        plot.subtitle = element_text(
          face = "italic", size = FIG_SUBTITLE_SIZE, colour = "grey30"
        )
      )
    )
  save_panel(panel, file.path(dirs$rpt, "panel"), width = 230, height = 190)

  write_leaf_workbook(
    file.path(dirs$dat, "source_data.xlsx"),
    list(
      contrasts_proteins = prot,
      contrasts_pathways = contrasts_by_space$pathways,
      contrasts_modules = contrasts_by_space$modules,
      variance_partition = vp
    )
  )
  invisible(contrasts_by_space)
}

# Global ordination: PCA with subject trajectories plus a PERMANOVA respecting
# repeated measures via strata = subject.
run_global_ordination_leaf <- function(bundle, out_dir) {
  pacman::p_load(vegan)
  dirs <- leaf_dirs(out_dir)
  da <- readRDS(assoc_paths()$norm)
  source(here("functions", "shared_pca.R"))

  pca <- run_pca(bundle$imputed_proteins, da$metadata, log_transform = FALSE)

  d <- t(bundle$imputed_proteins)
  md <- da$metadata[match(rownames(d), da$metadata$Col_ID), ]
  set.seed(1)
  perm <- vegan::adonis2(
    d ~ Group * Timepoint,
    data = md, method = "euclidean",
    strata = md$Subject_ID, permutations = 999, by = "terms"
  )
  perm_tab <- tibble::rownames_to_column(as.data.frame(perm), "term")

  p_pca <- build_pca_panel(pca$scores, pca$var_exp)
  p_tab <- gridExtra::tableGrob(
    perm_tab |>
      dplyr::filter(term %in% c("Group", "Timepoint", "Group:Timepoint")) |>
      dplyr::transmute(term,
        R2 = round(R2, 3), F = round(.data$F, 2),
        p = signif(.data$`Pr(>F)`, 2)
      ),
    rows = NULL, theme = gridExtra::ttheme_minimal(base_size = 8)
  )

  panel <- p_pca + patchwork::wrap_elements(p_tab) +
    patchwork::plot_annotation(
      title = "Global ordination and PERMANOVA",
      subtitle = "PERMANOVA on all samples, strata = subject. Descriptive.",
      theme = theme(
        plot.title = element_text(face = "bold", size = FIG_TITLE_SIZE),
        plot.subtitle = element_text(
          face = "italic", size = FIG_SUBTITLE_SIZE, colour = "grey30"
        )
      )
    )
  save_panel(panel, file.path(dirs$rpt, "panel"), width = 230, height = 120)

  write_leaf_workbook(
    file.path(dirs$dat, "source_data.xlsx"),
    list(pca_scores = pca$scores, permanova = perm_tab)
  )
  invisible(perm_tab)
}

# Group HR-vs-LR leaf: one contrast, one method, three feature-space volcanoes.
# limma is the leaner per-contrast snapshot (with a 03_DEP reconciliation);
# mixed_slice reads the cached global fit.
run_group_leaf <- function(bundle, contrast, method, out_dir) {
  dirs <- leaf_dirs(out_dir)
  group_of <- assoc_group_of(bundle$meta)

  res_by_space <- if (method == "limma") {
    lapply(bundle$feature_sets, function(mat) {
      group_limma_contrast(mat, bundle$meta, contrast, group_of)
    })
  } else {
    cache <- readRDS(global_cache_path())
    lapply(cache, function(tab) {
      tab |>
        filter(.data$contrast == !!contrast) |>
        transmute(feature, effect = estimate, t, p, bh)
    })
  }

  sheets <- res_by_space
  if (method == "limma") {
    recon <- res_by_space$proteins |>
      inner_join(dep_reference(contrast), by = "feature")
    sheets$reconciliation_vs_03DEP <- recon
  }

  write_leaf_workbook(file.path(dirs$dat, "source_data.xlsx"), sheets)
  panel <- build_group_leaf_panel(
    res_by_space, bundle$gene_map, contrast, method
  )
  save_panel(panel, file.path(dirs$rpt, "panel"), width = 230, height = 95)
  invisible(res_by_space)
}

# Phenotype leaf: feature ~ continuous adaptation, one phase, six traits, three
# feature spaces. Volcano leads with the composite hypertrophy trait.
run_phenotype_leaf <- function(bundle, phase, out_dir) {
  dirs <- leaf_dirs(out_dir)
  res_by_space <- lapply(bundle$feature_sets, function(mat) {
    associate_limma(mat, bundle$meta_da, bundle$pheno, PHENO_TRAITS,
      phases = phase
    )
  })

  vs <- lapply(names(res_by_space), function(sp) {
    lead <- res_by_space[[sp]] |>
      filter(trait == "comp_hypertrophy") |>
      transmute(feature, effect = beta, t, p, bh)
    assoc_volcano(lead,
      gene_map = if (sp == "proteins") bundle$gene_map else NULL,
      title = SPACE_LABELS_A[[sp]], effect_lab = "slope vs hypertrophy"
    )
  })
  cap <- if (phase == "acute") {
    "Acute signal reflects muscle-damage repair, not hypertrophy (Damas 2016)."
  } else {
    "In-sample association; BH within phase x trait."
  }
  panel <- wrap_plots(vs, nrow = 1) +
    plot_annotation(
      title = sprintf("Phenotype association — %s phase", phase),
      subtitle = cap,
      theme = theme(
        plot.title = element_text(face = "bold", size = FIG_TITLE_SIZE),
        plot.subtitle = element_text(
          face = "italic", size = FIG_SUBTITLE_SIZE, colour = "grey30"
        )
      )
    )
  save_panel(panel, file.path(dirs$rpt, "panel"), width = 230, height = 95)

  write_leaf_workbook(
    file.path(dirs$dat, "source_data.xlsx"),
    setNames(res_by_space, paste0("assoc_", names(res_by_space)))
  )
  invisible(res_by_space)
}

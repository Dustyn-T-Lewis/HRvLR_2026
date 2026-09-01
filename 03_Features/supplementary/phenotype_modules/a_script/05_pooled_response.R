#!/usr/bin/env Rscript
# The pooled all-subject response, and the pathway layer over both it
# and the phenotype associations.
#
# The nine canonical contrasts split every time effect by arm; nothing
# ever tested the average response across all 16 subjects. Two pooled
# contrasts on the same means model, blocking and consensus correlation
# as the canonical fit: Training_All, the mean of the two T2 cells minus
# the mean of the two T1 cells, and Acute_All, T3 minus T2 likewise.
#
# The association analyses upstream need no such contrast (they use
# per-subject deltas); this is the backdrop that says what training and
# an acute bout do on average, with double the n of the within-arm
# splits. Pathways run through both engines: limma::fry, whose rotation
# null carries inter-gene correlation and is the test this repository
# trusts, and preranked fgsea for effect direction. fgsea also scores
# the phenotype associations — proteins ranked by their delta-phenotype
# Spearman rho — where its gene-permutation null is the only one
# available; those rows are descriptive and labelled so.
#
# Read 05_pooled_dep.csv, 05_pooled_pathways.csv (fry + fgsea per pooled
# contrast) and 05_fgsea_pheno.csv (config x phenotype NES). A null for
# the phenotype rows looks like no set at padj < 0.05; the pooled
# contrasts themselves are expected to be strongly non-null, because a
# main effect of exercise is not in question.

pacman::p_load(here, dplyr, tidyr, readr, tibble, purrr, limma, fgsea)
source(here("functions", "feature_contrasts.R"))
source(here("functions", "shared_pathway_utils.R"))
source(here("functions", "shared_hlm.R"))
source(here("functions", "blood_index_model.R"))
source(here(
  "03_Features", "supplementary", "galamm_pilot", "a_script", "pilot_helpers.R"
))
source(here(
  "03_Features", "supplementary", "phenotype_modules", "a_script",
  "trait_helpers.R"
))

set.seed(42)
OUT <- here("03_Features", "supplementary", "phenotype_modules", "c_data")

POOLED <- c(
  Training_All = "(HR_T2 + LR_T2)/2 - (HR_T1 + LR_T1)/2",
  Acute_All = "(HR_T3 + LR_T3)/2 - (HR_T2 + LR_T2)/2"
)
PHENOS <- c(
  "comp_hypertrophy", "d_fcsa_I", "d_fcsa_II", "d_mcsa",
  "d_1rm_legpress", "d_1rm_ext"
)

mat <- protein_matrix()
meta <- feature_metadata()
parts <- feature_design(mat, meta)
cm <- limma::makeContrasts(contrasts = unname(POOLED), levels = parts$design)
colnames(cm) <- names(POOLED)

fit <- limma::lmFit(mat, parts$design,
  block = parts$block, correlation = parts$correlation
)
fit2 <- limma::eBayes(limma::contrasts.fit(fit, cm))

anno <- readRDS(
  here("02_Normalization", "c_data", "DAList_normalized.rds")
)$annotation

pooled <- map_dfr(colnames(cm), function(ct) {
  limma::topTable(fit2,
    coef = ct, number = Inf, adjust.method = "BH",
    sort.by = "none"
  ) |>
    rownames_to_column("feature") |>
    transmute(
      contrast = ct, feature = .data$feature,
      gene = anno$gene[match(.data$feature, anno$uniprot_id)],
      logFC = .data$logFC, t = .data$t, p = .data$P.Value,
      bh = .data$adj.P.Val
    )
})
write_csv(pooled, file.path(OUT, "05_pooled_dep.csv"))
message(paste(
  capture.output(print(as.data.frame(
    summarise(pooled,
      bh05 = sum(.data$bh < 0.05), .by = "contrast"
    )
  ))),
  collapse = "\n"
))

hallmark <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
sets_symbol <- split(hallmark$gene_symbol, hallmark$gs_name)
sets_uniprot <- lapply(sets_symbol, function(g) {
  anno$uniprot_id[anno$gene %in% g]
})

mat_complete <- mat[stats::complete.cases(mat), , drop = FALSE]
fry_res <- pathway_fry(
  mat_complete, sets_uniprot, parts$design, cm,
  block = parts$block, correlation = parts$correlation
)

rank_of <- function(stat, gene) {
  ok <- !is.na(gene) & !is.na(stat)
  tapply(stat[ok], gene[ok], function(x) x[which.max(abs(x))])
}
fgsea_pooled <- map_dfr(colnames(cm), function(ct) {
  rows <- filter(pooled, .data$contrast == ct)
  run_fgsea(sort(rank_of(rows$t, rows$gene)), sets_symbol) |>
    mutate(contrast = ct)
})

write_csv(
  full_join(
    rename(fry_res, fry_p = "p", fry_fdr = "fdr"),
    select(
      fgsea_pooled, "contrast", "pathway",
      fgsea_nes = "NES", fgsea_padj = "padj"
    ),
    by = c("contrast", "pathway")
  ),
  file.path(OUT, "05_pooled_pathways.csv")
)
message(paste(
  capture.output(print(as.data.frame(
    fry_res |>
      filter(.data$fdr < 0.05) |>
      count(.data$contrast, .data$direction)
  ))),
  collapse = "\n"
))

inp <- pilot_data()
prot_ids <- tibble(
  subject = sub("^(HR|LR)_", "", as.character(inp$meta$subject)),
  timepoint = as.character(inp$meta$timepoint)
)
pheno <- read_csv(
  here("00_input", "c_data", "phenotype.csv"),
  show_col_types = FALSE
) |>
  mutate(subject = sub("^(HR|LR)_", "", .data$subject))
gene_931 <- anno$gene[match(rownames(inp$mat), anno$uniprot_id)]

delta_for <- function(config) {
  have <- split(prot_ids$subject, prot_ids$timepoint)
  subjects <- intersect(have[[config[1]]], have[[config[2]]])
  key <- paste(prot_ids$subject, prot_ids$timepoint, sep = "_")
  out <- inp$mat[, match(paste(subjects, config[2], sep = "_"), key)] -
    inp$mat[, match(paste(subjects, config[1], sep = "_"), key)]
  colnames(out) <- subjects
  out
}

fgsea_pheno <- imap_dfr(
  list(training = c("T1", "T2"), acute = c("T2", "T3")),
  function(config, config_name) {
    dm <- delta_for(config)
    ph <- pheno[match(colnames(dm), pheno$subject), ]
    map_dfr(PHENOS, function(p) {
      rho <- suppressWarnings(cor(t(dm), ph[[p]],
        method = "spearman", use = "pairwise.complete.obs"
      ))[, 1]
      run_fgsea(sort(rank_of(rho, gene_931)), sets_symbol) |>
        mutate(config = config_name, phenotype = p)
    })
  }
)
write_csv(
  select(fgsea_pheno, -"leadingEdge"),
  file.path(OUT, "05_fgsea_pheno.csv")
)
message(paste(
  capture.output(print(as.data.frame(
    fgsea_pheno |>
      filter(.data$padj < 0.05) |>
      count(.data$config, .data$phenotype)
  ))),
  collapse = "\n"
))

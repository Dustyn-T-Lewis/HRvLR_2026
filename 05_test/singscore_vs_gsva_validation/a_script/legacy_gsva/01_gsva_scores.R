# Per-sample pathway activity via GSVA (Gaussian kcdf for log-intensity data).
if (!exists("expr")) source(here::here("05_test", "a_script", "setup.R"))

gsva_par <- GSVA::gsvaParam(expr, GENE_SETS, minSize = 10, maxSize = 500)
scores <- GSVA::gsva(gsva_par, verbose = FALSE)

write_csv(
  as.data.frame(scores) |> tibble::rownames_to_column("pathway"),
  file.path(DAT, "gsva_scores.csv")
)
cat(sprintf("GSVA: %d pathways x %d samples\n", nrow(scores), ncol(scores)))

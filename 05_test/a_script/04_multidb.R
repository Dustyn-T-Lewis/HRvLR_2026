# Score GSVA and ssGSEA on all five databases, fit the group x timepoint mixed
# model per pathway per method, and correct WITHIN each database (each database
# has its own size and redundancy, so pooled BH is the wrong reference). Keeps
# nominal p alongside the within-database FDR.
if (!exists("expr")) source(here::here("05_test", "a_script", "setup.R"))
pacman::p_load(GSVA, lme4, lmerTest)

pw_full <- build_pathway_collection(
  min_size = 15, max_size = 500, include_goslim = TRUE, exclude_variants = TRUE
)

getm <- function(eff, g, t) {
  col <- paste0(g, "_", t)
  if (col %in% names(eff)) eff[[col]] else NA_real_
}

lmm_dir <- function(scores) {
  long <- scores |>
    as.data.frame() |>
    tibble::rownames_to_column("pathway") |>
    pivot_longer(-pathway, names_to = "sample", values_to = "score") |>
    left_join(meta, by = "sample")
  bind_rows(lapply(rownames(scores), function(pw) {
    d <- long[long$pathway == pw, ]
    tryCatch(
      {
        a <- anova(lmerTest::lmer(score ~ group * timepoint + (1 | subject), data = d))
        eff <- d |>
          group_by(group, timepoint) |>
          summarise(mu = mean(score), .groups = "drop") |>
          pivot_wider(names_from = c(group, timepoint), values_from = mu)
        slope <- (getm(eff, "HR", "T3") - getm(eff, "HR", "T1")) -
          (getm(eff, "LR", "T3") - getm(eff, "LR", "T1"))
        tibble(pathway = pw, p_interaction = a["group:timepoint", "Pr(>F)"], slope = slope)
      },
      error = function(e) tibble(pathway = pw, p_interaction = NA_real_, slope = NA_real_)
    )
  }))
}

s_gsva <- gsva(gsvaParam(expr, pw_full, minSize = 10, maxSize = 500), verbose = FALSE)
s_ssgsea <- gsva(ssgseaParam(expr, pw_full, minSize = 10, maxSize = 500), verbose = FALSE)

res <- bind_rows(
  lmm_dir(s_gsva) |> mutate(method = "GSVA"),
  lmm_dir(s_ssgsea) |> mutate(method = "ssGSEA")
) |>
  mutate(database = classify_database(pathway)) |>
  group_by(method, database) |>
  mutate(fdr_within = p.adjust(p_interaction, "BH")) |>
  ungroup()
write_csv(res, file.path(DAT, "multidb_results.csv"))

summ <- res |>
  group_by(method, database) |>
  summarise(
    n = sum(!is.na(p_interaction)),
    nominal = sum(p_interaction < 0.05, na.rm = TRUE),
    fdr = sum(fdr_within < 0.05, na.rm = TRUE), .groups = "drop"
  )
print(summ)
cat(sprintf("multi-DB: %d method-pathway models\n", nrow(res)))

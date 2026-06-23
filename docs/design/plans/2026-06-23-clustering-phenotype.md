# HRvLR Clustering → Phenotype Implementation Plan

**Goal:** Group the pi-gated HRvLR proteome with three pure-R engines, summarise each group by an eigengene, and link eigengenes to phenotype with a mixed model — WGCNA inferential, the rest corroborating or quarantined.

**Architecture:** Each engine is an independent script reading shared inputs via one sourced helper, writing the same two artifacts (`membership.csv`, `eigengene.csv`) to its `c_data/`. A single `d_integration` step consumes those for the LMM, permutation null, BH, cross-engine concordance, and ORA. Spec: `docs/design/specs/2026-06-23-hrvlr-clustering-phenotype-design.md`.

**Tech Stack:** R only, no Python. `here`, `pacman`; WGCNA, Mfuzz, mixOmics, lmerTest, clusterProfiler — each used through its own maintained API.

## Global Constraints

- Effective n = 16 subjects (HR 8 / LR 8), unbalanced — two HR subjects partial (S28: T1 only; S29: T2/T3 only); 45 samples are pseudo-replicated. Hypothesis-generating, not confirmatory.
- HR-only attrition (S28, S29) — carry as a caveat on HR-specific results.
- Clustering input AND ORA background = proteins with pi-score < 0.05 in any group (the same set).
- Pure R: no reticulate/Python anywhere.
- Scripts lean and vignette-grounded: read each package's GitHub/vignette, cite it in a short header, call the package's own functions, minimum code to produce the two artifacts. No defensive bloat, no narrative comments, no banner separators.
- Paths via `here::here()`; packages via `pacman::p_load(...)`. Match `03_DEP` conventions.
- Artifact contract (every engine writes to its `c_data/`):
  - `membership.csv`: `protein_id, group_id, membership_weight`
  - `eigengene.csv`: `sample_id, subject, group_arm, timepoint, group_id, ME`
- Phenotype resolved: per-subject, lead `COMP.HYPERTROPHY` + ΔT1→T2 components (fCSA-I/II, MyoVision fCSA-I, mCSA, 1-RM leg-press/ext) from `00_input/HRvLR_meta.csv`, join on `Subject_ID`. **No `group_arm` in models** (collinear with the group-defining composite). Timepoints: T1→T2 chronic, T3 acute. Selection caveat: pi-gate includes HR-vs-LR contrasts → permutation null is the guard, layer stays hypothesis-generating.

---

### Task 1: Shared input helper

**Files:**
- Create: `05_Clustering/_shared_inputs.R`
- Reads: `03_DEP/c_data/01_limma_DAList.rds`, `03_DEP/c_data/03_combined_results.csv`

**Interfaces:**
- Produces: `load_clustering_inputs()` returning a list with
  - `abund`: matrix, pi-gated proteins × 45 samples (imputed, normalized abundances)
  - `meta`: data.frame `sample_id, subject, group_arm, timepoint`
  - `gap`: matrix, pi-gated proteins × 3 (HR−LR limma logFC at T1/T2/T3)
  - `pi_set`: character vector of pi-gated protein ids (the ORA background)

- [ ] **Step 1: Inspect the source objects** — confirm where pi-scores and per-contrast logFCs live.

Run: `Rscript -e 'x <- readRDS(here::here("03_DEP/c_data/01_limma_DAList.rds")); str(x, max.level = 2)'`
Expected: see the abundance matrix, sample metadata, and per-contrast stats slots; note the exact names.

- [ ] **Step 2: Write `load_clustering_inputs()`** — pi-gate (pi-score < 0.05 in any group), subset the abundance matrix to that set, assemble `meta` from the DAList sample annotation, and build `gap` from the three HR−LR contrast logFCs. One function, returns the list above. Header cites that it reads `03_DEP` outputs.

- [ ] **Step 3: Verify shapes**

Run: `Rscript -e 'source(here::here("05_Clustering/_shared_inputs.R")); i <- load_clustering_inputs(); cat(nrow(i$abund), ncol(i$abund), nrow(i$gap), ncol(i$gap), length(i$pi_set), "\n"); stopifnot(ncol(i$abund) == 45, ncol(i$gap) == 3, nrow(i$abund) == length(i$pi_set), !anyNA(i$abund))'`
Expected: prints `<N> 45 <N> 3 <N>`; no error.

- [ ] **Step 4: Commit**

```bash
git add 05_Clustering/_shared_inputs.R
git commit -m "add shared input loader for clustering stage"
```

---

### Task 2: WGCNA paired (PRIMARY)

**Files:**
- Create: `05_Clustering/a_wgcna_paired/a_script/01_wgcna_paired.R`
- Write: `a_wgcna_paired/c_data/{membership.csv,eigengene.csv}`, `a_wgcna_paired/b_reports/{soft_threshold.pdf,dendro_modules.pdf}`

**Interfaces:**
- Consumes: `load_clustering_inputs()` (`abund`, `meta`)
- Produces: the two-artifact contract; `group_id` = module color/label.

**Vignette to read first:** WGCNA tutorials (horvath lab; github.com/cran/WGCNA) and the paired-design modification (Li et al. 2018, *Sci Rep*). Canonical: within-subject center `abund` → `pickSoftThreshold()` → `blockwiseModules()` (signed) → `moduleEigengenes()`.

- [ ] **Step 1: Write the script** — subtract each subject's across-timepoint mean from `abund` (within-subject centering, the step that stops modules encoding subject identity), transpose to samples × proteins. `pickSoftThreshold` (pick at scale-free R²>0.8), `blockwiseModules(networkType="signed", minModuleSize=20, mergeCutHeight=0.25, deepSplit=2)`. Eigengene = `moduleEigengenes()$eigengenes`, sign-fixed so positive ME = higher mean of positive-loading proteins, reshaped to the contract; membership = module assignment (weight 1). Report module-count sensitivity to power and `mergeCutHeight`.

- [ ] **Step 2: Run it**

Run: `Rscript 05_Clustering/a_wgcna_paired/a_script/01_wgcna_paired.R`
Expected: writes both CSVs + two PDFs; prints chosen power and module count.

- [ ] **Step 3: Verify artifacts**

Run: `Rscript -e 'e <- read.csv(here::here("05_Clustering/a_wgcna_paired/c_data/eigengene.csv")); m <- read.csv(here::here("05_Clustering/a_wgcna_paired/c_data/membership.csv")); stopifnot(all(c("sample_id","subject","group_arm","timepoint","group_id","ME") %in% names(e)), nrow(e) >= 45, all(c("protein_id","group_id","membership_weight") %in% names(m)))'`
Expected: no error.

- [ ] **Step 4: Commit**

```bash
git add 05_Clustering/a_wgcna_paired/
git commit -m "add within-subject-centered paired wgcna as primary clustering engine"
```

---

### Task 3: Mfuzz on HR−LR gap (companion — temporal shape)

**Files:**
- Create: `05_Clustering/b_mfuzz_gap/a_script/01_mfuzz_gap.R`
- Write: `b_mfuzz_gap/c_data/{membership.csv,eigengene.csv}`, `b_mfuzz_gap/b_reports/{gap_clusters.pdf,shape_profiles.pdf}`

**Interfaces:**
- Consumes: `load_clustering_inputs()` (uses `gap` and `abund`/`meta`)
- Produces: the two-artifact contract; `group_id` = Mfuzz cluster.

**Engine:** `e1071::cmeans` — the fuzzy-c-means Mfuzz wraps. `Mfuzz` itself imports `tcltk`/`tkWidgets` and won't load without XQuartz (absent here), so use `cmeans` directly and reimplement the thin Mfuzz helpers: `standardise` = z-score per protein (`scale` on the transpose); `mestimate` = the Schwämmle & Jensen (2010) fuzzifier formula `m = 1 + (1418/n + 22.05) * D^-2 + (12.33/n + 0.243) * D^(-0.0406*log(n) - 0.1134)` where n = #proteins, D = #dims (3); `Dmin` = min inter-centroid distance across candidate `c`; `acore` = membership thresholding. Read the Mfuzz vignette/source (bioconductor.org/packages/Mfuzz) to match these.

- [ ] **Step 1: Write the script** — cluster the 3-point `gap` trajectory. `m` from `mestimate`; `c` chosen via `Dmin` capped at the geometry (~4–6); membership cutoff ≥ 0.5; loop ≥25 seeds and report assignment stability. Eigengene = first PC of each cluster's per-sample abundances (sign-fixed so positive ME = higher mean of positive-loading proteins). Membership from `acore`. Also emit the descriptive group-mean 6-point shape plot (illustration only). Header notes gap = inference-eligible, shape = illustration.

- [ ] **Step 2: Run it**

Run: `Rscript 05_Clustering/b_mfuzz_gap/a_script/01_mfuzz_gap.R`
Expected: writes both CSVs + two PDFs; prints chosen `c`, `m`, seed-stability.

- [ ] **Step 3: Verify artifacts**

Run: `Rscript -e 'e <- read.csv(here::here("05_Clustering/b_mfuzz_gap/c_data/eigengene.csv")); stopifnot(all(c("sample_id","subject","group_arm","timepoint","group_id","ME") %in% names(e)), !anyNA(e$ME))'`
Expected: no error.

- [ ] **Step 4: Commit**

```bash
git add 05_Clustering/b_mfuzz_gap/
git commit -m "add mfuzz gap-trajectory clustering as temporal companion"
```

---

### Task 4: Integration — concordance + ORA (phenotype-free, runs now)

**Files:**
- Create: `05_Clustering/d_integration/a_script/01_concordance_ora.R`
- Write: `d_integration/c_data/{concordance.csv,ora_results.csv}`, `d_integration/b_reports/{concordance_heatmap.pdf,ora_dotplots.pdf}`

**Interfaces:**
- Consumes: each engine's `membership.csv` (WGCNA, Mfuzz); `load_clustering_inputs()$pi_set` as background
- Produces: cross-engine membership agreement; ORA tables per group.

**Vignette to read first:** clusterProfiler book (github.com/YuLab-SMU/clusterProfiler). Use `enricher()` with `universe = pi_set`; gene sets GO BP/CC/MF + Reactome; BH within each source; **no `simplify()`**.

- [ ] **Step 1: Write the script** — load each available engine's membership; compute pairwise cross-engine concordance (adjusted Rand index + Jaccard on group overlaps) → `concordance.csv` + heatmap. Run ORA per group of every engine against `pi_set` background → `ora_results.csv` + dotplots.

- [ ] **Step 2: Run it**

Run: `Rscript 05_Clustering/d_integration/a_script/01_concordance_ora.R`
Expected: writes the two CSVs + two PDFs.

- [ ] **Step 3: Verify**

Run: `Rscript -e 'c <- read.csv(here::here("05_Clustering/d_integration/c_data/concordance.csv")); stopifnot(nrow(c) > 0); o <- read.csv(here::here("05_Clustering/d_integration/c_data/ora_results.csv")); stopifnot("p.adjust" %in% names(o))'`
Expected: no error.

- [ ] **Step 4: Commit**

```bash
git add 05_Clustering/d_integration/a_script/01_concordance_ora.R 05_Clustering/d_integration/c_data 05_Clustering/d_integration/b_reports
git commit -m "add cross-engine concordance and ora with pi-gated background"
```

---

### Task 5a: Phenotype table (per subject)

**Files:**
- Create: `05_Clustering/d_integration/a_script/02_phenotype_table.R`
- Reads: `00_input/HRvLR_meta.csv`
- Write: `d_integration/c_data/phenotype.csv`

**Interfaces:**
- Produces: per-subject phenotype, one row per subject (`subject`, `group_arm`, one
  column per trait), consumed by Task 5b.

- [ ] **Step 1** — read the long meta; compute per-subject ΔT1→T2 for `fCSA_Type_I_Pre`, `fCSA_Type_II_Pre`, `MyoVision_fCSA_Type_I__Pre`, `mCSA_Pre`, `X1RM_Leg_Pre`, `X1RM._Ext_Pre` (T2 value − T1 value); take `COMP.HYPERTROPHY` from the T2 row (`readr::parse_number`) as `comp_hypertrophy` (the lead). `subject` = `Subject_ID`, matching the proteome. **No `group_arm` enters the models** but keep it as a label column. Subjects missing a needed timepoint (S28: T1 only → Δ all NA; S29: T2/T3 only → Δ all NA, comp present) get NA for those traits; do not impute. Report per-trait non-missing n.

- [ ] **Step 2: Run** `Rscript 05_Clustering/d_integration/a_script/02_phenotype_table.R`. Expected: writes `phenotype.csv`, prints per-trait n.

- [ ] **Step 3: Verify**
Run: `Rscript -e 'p <- read.csv(here::here("05_Clustering/d_integration/c_data/phenotype.csv")); stopifnot(nrow(p)==16, "comp_hypertrophy" %in% names(p), "subject" %in% names(p)); cat("traits:", ncol(p)-2, " complete comp:", sum(!is.na(p$comp_hypertrophy)), "\n")'`
Expected: 16 rows; prints trait count and ~15 complete comp.

- [ ] **Step 4: Commit** `git add 05_Clustering/d_integration/a_script/02_phenotype_table.R 05_Clustering/d_integration/c_data/phenotype.csv && git commit -m "build per-subject phenotype table for clustering linkage"`

---

### Task 5b: Three phenotype models + permutation null

**Files:**
- Create: `05_Clustering/d_integration/a_script/03_phenotype_models.R`
- Write: `d_integration/c_data/{phenotype_models.csv,perm_null.csv}`, `d_integration/b_reports/{pheno_grid.pdf,lead_forest.pdf}`

**Interfaces:**
- Consumes: `phenotype.csv` (5a); each engine's `eigengene.csv`.
- Produces: the trait × module effect grid for three models + the falsification null.

**Vignette to read first:** lmerTest; a partial-R² helper (`r2glmm`/`partR2` if installed, else report marginal R²).

- [ ] **Step 1** — for each engine (WGCNA primary, Mfuzz secondary) × module × trait, fit three models on complete cases (no `group_arm`):
  1. *tracks*: `lmerTest::lmer(ME ~ pheno + timepoint + (1|subject))` over all samples → standardised β, partial/marginal R², p.
  2. *baseline*: `lm(pheno ~ T1_ME)` (per-subject T1 eigengene) → β, R², p.
  3. *acute*: `lm(pheno ~ (T3_ME − T2_ME))` (per-subject acute shift) → β, R², p.
  Collect to `phenotype_models.csv` (`engine, module, trait, model, n, beta_std, r2, p, p_bh`); **BH across the whole grid**. Skip cells with n < 6 (report skipped).
- [ ] **Step 2** — permutation null: permute the phenotype vector across subjects (B=1000), recompute the max |stat| per model; store the null in `perm_null.csv` and a `p_perm` for the observed max. Subjects with NA phenotype are excluded before permuting (permute only the observed values).
- [ ] **Step 3** — plots: `pheno_grid.pdf` = trait × module effect heatmaps (one facet per model, signed β, BH-significant cells marked), `lead_forest.pdf` = forest of the WGCNA-module effects for `comp_hypertrophy` (the lead) with the permutation-null band. State the selection caveat + hypothesis-generating label on the figures.
- [ ] **Step 4: Run** `Rscript 05_Clustering/d_integration/a_script/03_phenotype_models.R`. Expected: both CSVs + both PDFs; prints whether any cell survives BH and whether the observed max beats the null.
- [ ] **Step 5: Verify**
Run: `Rscript -e 'm <- read.csv(here::here("05_Clustering/d_integration/c_data/phenotype_models.csv")); stopifnot(all(c("engine","module","trait","model","beta_std","p_bh") %in% names(m)), all(c("tracks","baseline","acute") %in% m$model))'`
Expected: no error.
- [ ] **Step 6: Commit** `git add 05_Clustering/d_integration/a_script/03_phenotype_models.R 05_Clustering/d_integration/c_data 05_Clustering/d_integration/b_reports && git commit -m "add three module-phenotype models with subject-permutation null"`

---

### Task 6: Supervised sPLS/DIABLO (BLOCKED on phenotype; quarantine)

**Precondition:** phenotype inputs confirmed. Labeled exploratory throughout — never the headline.

**Files:**
- Create: `05_Clustering/c_supervised/a_script/01_spls_diablo.R`
- Write: `c_supervised/c_data/{membership.csv,cv_error.csv}`, `c_supervised/b_reports/spls_components.pdf`

**Interfaces:**
- Consumes: `load_clustering_inputs()$abund`; phenotype table
- Produces: candidate proteins (membership contract) + honest nested-CV error.

**Vignette to read first:** mixOmics vignettes (mixomics.org; github.com/mixOmicsTeam/mixOmics). `spls()`/`block.splsda()`, `tune.spls()`, `perf()` with nested CV.

- [ ] **Step 1: Write the script** — sPLS (proteome vs phenotype) or DIABLO if multi-block; nested CV via `tune.spls`/`perf`; report CV error honestly in `cv_error.csv`. Selected proteins → `membership.csv`. Header and report state: exploratory, n=15 overfits, not confirmatory.

- [ ] **Step 2: Run it**

Run: `Rscript 05_Clustering/c_supervised/a_script/01_spls_diablo.R`
Expected: writes both CSVs + the PDF; prints CV error.

- [ ] **Step 3: Verify**

Run: `Rscript -e 'm <- read.csv(here::here("05_Clustering/c_supervised/c_data/membership.csv")); stopifnot(all(c("protein_id","group_id","membership_weight") %in% names(m)))'`
Expected: no error.

- [ ] **Step 4: Commit**

```bash
git add 05_Clustering/c_supervised/
git commit -m "add exploratory supervised spls with nested cv, quarantined"
```

---

## Self-Review

**Spec coverage:** input/background → Task 1, 4; WGCNA primary + centering → Task 2; Mfuzz gap + shape → Task 3; concordance + ORA → Task 4; LMM + permutation null + BH → Task 5; supervised quarantine → Task 6; two open inputs gate 5/6 explicitly. MEFISTO removed (pure-R constraint). All spec sections covered.

**Placeholder scan:** no "TBD/handle edge cases"; the two blocked tasks have explicit preconditions, not vague gaps. Engine cores intentionally specify canonical functions + params + vignette rather than fabricated code, per the vignette-first constraint.

**Type consistency:** the two-artifact column contract is identical across Tasks 2, 3, 6 and is what 4/5 consume; `group_id` is the per-engine grouping key throughout; `pi_set` is both clustering input and ORA background.

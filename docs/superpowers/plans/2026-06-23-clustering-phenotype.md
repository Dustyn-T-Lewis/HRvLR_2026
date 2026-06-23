# HRvLR Clustering → Phenotype Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Group the pi-gated HRvLR proteome with three pure-R engines, summarise each group by an eigengene, and link eigengenes to phenotype with a mixed model — WGCNA inferential, the rest corroborating or quarantined.

**Architecture:** Each engine is an independent script reading shared inputs via one sourced helper, writing the same two artifacts (`membership.csv`, `eigengene.csv`) to its `c_data/`. A single `d_integration` step consumes those for the LMM, permutation null, BH, cross-engine concordance, and ORA. Spec: `docs/superpowers/specs/2026-06-23-hrvlr-clustering-phenotype-design.md`.

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
- Phenotype-blind engines (Tasks 1–4) run now. Supervised (Task 6) and the LMM half of integration (Task 5) are blocked on two confirmed inputs: (a) phenotype variable(s) + per-timepoint vs per-subject, (b) timepoint coding / T3 handling (default: T1→T2 training axis, T3 a separate acute level).

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

### Task 5: Integration — LMM + permutation null (BLOCKED on phenotype inputs)

**Precondition:** phenotype variable(s) confirmed + per-timepoint vs per-subject decided; timepoint/T3 coding confirmed. Do not start until both are set.

**Files:**
- Create: `05_Clustering/d_integration/a_script/02_phenotype_lmm.R`
- Write: `d_integration/c_data/{lmm_results.csv,perm_null.csv}`, `d_integration/b_reports/lmm_forest.pdf`

**Interfaces:**
- Consumes: every engine's `eigengene.csv`; the phenotype table
- Produces: per-group effect sizes + the falsification test.

**Vignette to read first:** lmerTest + a partial-R² package (`partR2` or `r2glmm`).

- [ ] **Step 1: Write the script** — for each module of the **primary** (WGCNA) engine: `lmer(ME ~ phenotype + group_arm + timepoint + (1|subject))` via lmerTest; standardised β + partial R². BH across the primary engine's modules only; companion engines reported as concordant/discordant, not added to the correction. Permutation null: permute phenotype across the 16 subjects (subject-level, keeping each subject's actual timepoint block intact — two HR subjects have <3), recompute the largest |group–phenotype| association, B = 1000 → `perm_null.csv`. Forest plot of effects with the null band.

- [ ] **Step 2: Run it**

Run: `Rscript 05_Clustering/d_integration/a_script/02_phenotype_lmm.R`
Expected: writes both CSVs + the forest PDF; prints whether the top association beats the null.

- [ ] **Step 3: Verify**

Run: `Rscript -e 'l <- read.csv(here::here("05_Clustering/d_integration/c_data/lmm_results.csv")); stopifnot(all(c("group_id","beta_std","partial_r2","p_bh") %in% names(l)))'`
Expected: no error.

- [ ] **Step 4: Commit**

```bash
git add 05_Clustering/d_integration/a_script/02_phenotype_lmm.R 05_Clustering/d_integration/c_data 05_Clustering/d_integration/b_reports
git commit -m "add module-phenotype lmm with subject-permutation null"
```

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

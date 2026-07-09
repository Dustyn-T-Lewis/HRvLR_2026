# DEP stage correctness audit — A_HRvLR_2026

Read-only audit. Branch `refactor/hrvlr-cleanup`. No files modified. proteoDA
internals (`add_design`, `add_contrasts`, `fit_limma_model`,
`extract_DA_results`, `DAList`, `validate_DAList`) were inspected via the
installed package (`/Library/Frameworks/R.framework/Versions/4.6/Resources/library/proteoDA`)
to ground checks 2 and 3 in the actual executed logic, not assumption.

## Check 1 — Contrast algebra

File: `03_DEP/contrasts.R:4-14`.

Six cell means exist in the design (verified in Check 2): `HR_T1, HR_T2,
HR_T3, LR_T1, LR_T2, LR_T3`. Each of the nine contrasts parses as:

| # | Name | Left | Right | Resolves to |
|---|---|---|---|---|
| 1 | `Training_HR` | `HR_T2` | `- HR_T1` | `HR_T2 - HR_T1` |
| 2 | `Training_LR` | `LR_T2` | `- LR_T1` | `LR_T2 - LR_T1` |
| 3 | `Acute_HR` | `HR_T3` | `- HR_T2` | `HR_T3 - HR_T2` |
| 4 | `Acute_LR` | `LR_T3` | `- LR_T2` | `LR_T3 - LR_T2` |
| 5 | `Baseline_HRvLR` | `HR_T1` | `- LR_T1` | `HR_T1 - LR_T1` |
| 6 | `Trained_HRvLR` | `HR_T2` | `- LR_T2` | `HR_T2 - LR_T2` |
| 7 | `Acute_HRvLR` | `HR_T3` | `- LR_T3` | `HR_T3 - LR_T3` |
| 8 | `Training_Interaction` | `(HR_T2 - HR_T1)` | `- (LR_T2 - LR_T1)` | `HR_T2 - HR_T1 - LR_T2 + LR_T1` |
| 9 | `Acute_Interaction` | `(HR_T3 - HR_T2)` | `- (LR_T3 - LR_T2)` | `HR_T3 - HR_T2 - LR_T3 + LR_T2` |

All nine are linear combinations of only the six named cell means, each
coefficient is ±1, and each contrast's coefficients sum to zero (valid
contrast, not a cell-mean estimate). The two interaction terms expand
correctly to the difference-of-differences form required for a Group ×
Time interaction: contrast 8 tests whether the HR training response
(`T2−T1`) differs from the LR training response; contrast 9 does the same
for the acute bout (`T3−T2`). No sign errors — both interactions are
written `HR_delta − LR_delta`, consistently. Every name used
(`HR_T1…LR_T3`) is a design column (Check 2). `contrasts.R:16` sets
`PI_THRESH <- 0.05`, used only as a named constant, not part of the
algebra.

`03_DEP/a_non_imputed/a_script/01_run_dep.R:103` and
`03_DEP/b_imputed/a_script/01_run_dep_imputed.R:44` both call
`add_contrasts(dal, contrasts_vector = HRVLR_CONTRASTS)`, which internally
runs `limma::makeContrasts()` (proteoDA `add_contrasts`, inspected
directly) — this would itself throw a parse error at runtime if any
contrast failed to resolve against the design-matrix columns, and
`add_contrasts` additionally pre-checks `contrast_groups %in%
colnames(design_matrix)` before calling `makeContrasts`, aborting with a
named list of any invalid group. No such failure is possible given the
column names established in Check 2.

**Verdict: PASS**

## Check 2 — Design + block

Non-imputed arm, `03_DEP/a_non_imputed/a_script/01_run_dep.R:71-76,99,103,107`:
```r
meta$group <- factor(meta$group,
  levels = c("HR_T1", "HR_T2", "HR_T3", "LR_T1", "LR_T2", "LR_T3"))
...
dal <- add_design(dal, "~ 0 + group + (1 | subject)")
...
dal <- add_contrasts(dal, contrasts_vector = HRVLR_CONTRASTS)
...
dal <- fit_limma_model(dal)
```

Imputed arm, `03_DEP/b_imputed/a_script/01_run_dep_imputed.R:38-45`:
```r
dal$metadata$group <- factor(
  dal$metadata$Group_Time,
  levels = c("HR_T1", "HR_T2", "HR_T3", "LR_T1", "LR_T2", "LR_T3"))
dal$metadata$subject <- dal$metadata$Subject_ID
dal <- add_design(dal, "~ 0 + group + (1 | subject)")
dal <- add_contrasts(dal, contrasts_vector = HRVLR_CONTRASTS)
dal <- fit_limma_model(dal)
```

Both arms use the identical design-formula string, the identical
`levels =` ordering for `group`, and identical `subject` block column
name. Both source `HRVLR_CONTRASTS` from `03_DEP/contrasts.R` rather than
re-declaring or re-ordering it — no drift between arms.

Traced proteoDA internals to confirm the formula does what the comments
claim:

- `add_design()` (`proteoDA::add_design`, inspected): splits off the
  `(1 | subject)` term, fits `model.matrix(~0 + group, data = metadata)`,
  strips the `group` prefix from column names
  (`str_remove_all(colnames, "^group")`), and stores `random_factor =
  "subject"`. Because `group` is a factor with levels
  `c("HR_T1","HR_T2","HR_T3","LR_T1","LR_T2","LR_T3")`, the resulting
  design-matrix columns after prefix-stripping are exactly
  `HR_T1, HR_T2, HR_T3, LR_T1, LR_T2, LR_T3` in that order — matching
  every name used in `contrasts.R` (Check 1).
- `fit_limma_model()` (`proteoDA::fit_limma_model`, inspected): pulls
  `block <- DAList$metadata[, rf]` where `rf = design$random_factor =
  "subject"` (i.e., `Subject_ID`, per-subject not per-group-time), runs
  `limma::duplicateCorrelation(data, design_matrix, block = block)`,
  refits `lmFit(..., block = block, correlation = intra_block_cor)`,
  then — only if a contrast matrix is present —
  `contrasts.fit(fit, contrasts)`, and finally `eBayes(fit, robust =
  TRUE)`. eBayes is applied strictly after `contrasts.fit`, as required.
- `add_contrasts()` (`proteoDA::add_contrasts`, inspected): requires
  every contrast-referenced name to already be a design-matrix column
  name (verified above) before calling `limma::makeContrasts`.

Both scripts also guard the block itself:
`01_run_dep.R:78-80` — `if (any(is.na(meta$subject)) || any(meta$subject
== "")) stop(...)` — before fitting, so `duplicateCorrelation` can never
silently receive an `NA`/empty block level.

**Verdict: PASS.** No drift between arms; block is genuinely per-subject
(`Subject_ID`, e.g. `"HR_S03"`, unique within and across responder
groups — confirmed in `00_input/HRvLR_meta.csv`), not per-Group_Time.

## Check 3 — Naming and joins

Traced the metadata attachment chain from raw input through to the DEP
fit:

1. `00_input/HRvLR_meta.csv` carries `Col_ID, Subject_ID, Group,
   Timepoint, Group_Time` pre-computed (48 rows = 8 subjects × 2 groups ×
   3 timepoints), not regex-derived.
2. `01_Filtering/a_script/01_run_filtering.R:43-46`:
   `stopifnot("Sample mismatch" = setequal(colnames(intensity),
   metadata$Col_ID))`, then `intensity <- intensity[, metadata$Col_ID]`
   — reorders the matrix to metadata order explicitly (not a silent
   `match()` that could produce `NA` columns).
3. `DAList(data = int_mat, annotation = annot_df, metadata = metadata)`
   at `01_run_filtering.R:119` — proteoDA's `DAList()` constructor
   (inspected) calls `validate_DAList()`, which enforces
   `identical(colnames(x$data), rownames(x$metadata))` (exact order, not
   just set-equality) and aborts with an explicit hint
   ("same samples, different order — reorder metadata to match data
   columns") if violated. This is a hard runtime gate against the row/
   column-misalignment bug class noted in project memory
   (`feedback_dalist_misalignment.md`), not a place where a mismatch
   could pass silently.
4. `02_Normalization/a_script/01_run_normalization.R` only normalizes
   `dal$data` (cycloess) and writes `export_df <- bind_cols(annotation,
   as_tibble(dal$data))` — column order of `dal$data` is untouched, so
   `normalized.csv` columns keep filtering-stage order.
5. `03_DEP/a_non_imputed/a_script/01_run_dep.R:37-43,52-53,62-68`:
   `mat` is built from `normalized.csv` (`read_csv` preserves column
   order); `meta` is built from `dal_norm$metadata` (the *same*
   validated DAList's metadata, `readRDS` of
   `DAList_normalized.rds`), row order therefore identical to
   `dal$data` column order by construction (point 3's invariant is
   preserved through `saveRDS`/`readRDS`, which do not reorder).
   `rownames(meta_df) <- meta$sample_id` and the subsequent `DAList(data
   = mat, ..., metadata = meta_df)` call re-triggers the same
   `identical(colnames, rownames)` gate in `validate_DAList()` — so even
   if the two independent construction paths (CSV vs RDS) ever
   diverged in order, the script would hard-fail rather than silently
   join wrong subjects to wrong intensities.
   - Minor observation (not a defect): the script's own sanity check at
     line 83, `stopifnot(setequal(colnames(mat), meta$sample_id))`, is a
     set-equality check only, not an order check. It is not load-bearing
     for correctness because `DAList()`'s internal `validate_DAList()`
     performs the strict ordered check immediately afterward at line 90;
     flagging only for precision of intent (the comment/check implies
     less protection than actually exists).
6. Imputation scripts (`02_Normalization/imputation/a_script/{a_imp4p,
   b_mscoreutils,c_missforest,d_perseus}.R`) only replace `dal$data`
   values; `missForest`/`imp4p` transpose the matrix, sort/impute, then
   `dimnames(imp) <- dimnames(mat)` restores original row/column names
   and order before `stopifnot(identical(dim(imp), dim(mat)))`.
   `dal$metadata` is never touched by any of the four scripts, so the
   `01_run_dep_imputed.R` join described in Check 2 (`dal$metadata$group
   <-`, `dal$metadata$subject <-`) operates on the same validated,
   ordered metadata inherited from `DAList_normalized.rds`. No
   `match()`/join in the imputation scripts risks `NA` — the one
   `match()` present (`a_imp4p.R:17`,
   `cond <- factor(dal$metadata$Group_Time[match(colnames(mat),
   dal$metadata$Col_ID)])`) is redundant with the already-guaranteed
   alignment (used only to build the `imp4p` condition factor) but is
   not a source of drift since `colnames(mat)` and `dal$metadata$Col_ID`
   are, by the constructor invariant, the same set in the same order.
7. Factor levels: `time` is `factor(..., levels = c("T1","T2","T3"))`
   (`01_run_dep.R:70`) and `group` levels place all `*_T1` before
   `*_T2`/`*_T3` within each responder — T1 is consistently the
   pre/reference position everywhere it's declared. Because the model is
   cell-means (`~0 + group`, no intercept), level order does not change
   which coefficients are fit, only which design column appears first;
   contrast direction is instead pinned explicitly by the `A - B` syntax
   in `contrasts.R`, so T1-as-reference is enforced by the contrast
   strings, not by factor-level order — consistent with a cell-means
   design (no defect).

**Verdict: PASS.**

## Check 4 — Constants

- `PI_THRESH <- 0.05` (`03_DEP/contrasts.R:16`), documented in
  `01_run_dep.R:12-13` as `Pi = p^|logFC|`, threshold `Pi < 0.05`. Both
  DEP scripts compute `pi_score <- P.Value^abs(logFC)` identically
  (`01_run_dep.R:129`; `01_run_dep_imputed.R:53`) and both gate on
  `PI_THRESH`/`cfg$pi_thresh` sourced from the same constant — no
  re-declaration, no drift.
- BH FDR: `cfg$adj_method = "BH"` (`01_run_dep.R:26`) and
  `adj_method = "BH"` (`01_run_dep_imputed.R:46`) are both passed into
  `extract_DA_results()`, which (proteoDA source, inspected) forwards
  `adj_method` straight into `limma::topTable(..., adjust.method =
  adj_method)` to populate `adj.P.Val`. `pval_thresh = 0.10` is used
  consistently as the primary exploratory FDR cutoff in both arms, with
  `sig.FDR.05`/`sig.FDR.10` columns in the non-imputed summary
  (`01_run_dep.R:184-185,192-193`) as explicit fixed-cutoff companions —
  consistent with the FDR 0.10 rationale documented in
  `01_run_dep.R:113-116` and project memory.
- Fold-change direction: `sig_pi` is `+1` when `pi_score < thresh &
  logFC > 0` and `-1` when `logFC < 0`, in both scripts
  (`01_run_dep.R:130-133`; `01_run_dep_imputed.R:54-57`) — a positive
  `logFC` on e.g. `Training_HR = HR_T2 - HR_T1` means HR increased from
  T1 to T2, and this "left-minus-right, positive = increase in the
  left/first term" convention is applied uniformly across all nine
  contrasts (Check 1) and both arms. `lfc_thresh = 0` in both scripts
  means no additional fold-change magnitude filter is applied beyond
  the sign implied by `logFC`; this is a stated design choice
  (exploratory, FDR/Pi-gated) and is applied identically in both arms —
  not an inconsistency.

**Verdict: PASS.**

## Overall: PASS

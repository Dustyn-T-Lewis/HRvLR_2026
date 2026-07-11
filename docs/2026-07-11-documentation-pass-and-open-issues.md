# Documentation pass, and the code issues it did not fix

Date: 2026-07-11
Scope: the 16 Quarto narratives + the bibliography. **No analysis code was changed**, so
no committed number moved. Everything below that concerns *code* is open.

---

## What was fixed

### 1. The bibliography was fabricated. It is now Crossref-verified.

`03_DEP/a_non_imputed/a_script/references.bib` carried hallucinated author lists on four
entries and one outright wrong DOI, and it was **live** — `03_dep.qmd` cited `@smyth2005`
and `@xiao2014` into a rendered reference list.

| Entry | What was wrong |
|---|---|
| `xiao2014` | DOI `10.1093/bioinformatics/btu031` resolves to **InterProScan 5**, an unrelated paper. All 10 authors invented (incl. Kathryn Roeder, John Rinn — real scientists, uninvolved). Correct DOI: `10.1093/bioinformatics/btr671`, 7 authors, Xiao/Hsiao/Suresh/Chen/Wu/Wolf/Chen. |
| `smyth2005` | Pointed at Smyth **2004** (`10.2202/1544-6115.1027`), the *eBayes* paper — but was cited to justify `duplicateCorrelation`. Correct: Smyth, Michaud & Scott 2005, `10.1093/bioinformatics/bti270`. |
| `thurman_proteoDA_2023` | Correct DOI, entirely invented author list, wrong title and issue. |
| `arend_PRONE_2025` | Correct DOI, invented author list (fabricated Lukas Käll and Wolfgang Huber). |

**Action taken.** Every entry re-verified against Crossref. The bibliography now lives at
`shared/references.bib` (20 entries, all DOI-resolved), because stages 00–02 needed it and
burying it under `03_DEP` is why they had no citations at all. The old file is deleted.

### 2. Thirteen of 43 embedded figures were dead links

The panel-rename commits (`cf4fdea`, `143ccce`, `d393771`) broke the `![]()` paths in F02
(5), F04 (3), F05 (3) and F06 (2). Every F04/F05/F06 narrative rendered with holes. All 52
image links now resolve; ~8 previously-orphaned rendered panels are now embedded and
described (F04/F05 `panel_b_pattern_heatmap`, F06 `panel_b_nes_scatters`,
`supp_dendrogram`, `supp_soft_threshold`, the F03 top-30 supplements, the 05_test
prediction panels).

### 3. F03, F04 and F05 reported no results at all

They were pure method documents — no rho, no counts, no p-values, and therefore **no null
statement**, in exactly the figures where the result is null. They now state what was
found, pulled from the committed `c_data`.

### 4. Two narrative claims contradicted their own data

- **F06 said "leave-one-out Q² is not positive for the baseline adaptation associations."**
  It is positive in three cells, and the figure draws three white dots for them. The
  correct statement is that three of 84 cells have Q² > 0 — about what chance gives — and
  none survives BH (best q = 0.24).
- **F03 said its gate "is not the Xiao π-value (which is the *product*)".** It is. Taking
  log₁₀ of `p^|L| < 0.05` gives `|L| × (−log₁₀ p) > 1.301`. The exponential and product
  forms are the same rule; verified to select identical proteins on 200k random (p, logFC)
  pairs.

---

## Open issues — code, not documentation

These are now *disclosed* in the narratives where the reader meets them. None is *fixed*.
Fixing any of them changes published numbers and requires re-rendering.

### A. The π gate is not a significance filter, and it is stored in a column named `padj`

`01_run_dep.R:129` computes `pi_score = P.Value^abs(logFC)`, thresholded at 0.05.

- **489 π-calls against 12 proteins at BH < 0.05** — a 40× inflation.
- **71 of 489 (14.5%) have raw p ≥ 0.05.** Largest admitted raw p: **0.269**.
- The π set is therefore **not even a subset** of the nominally significant set.
- `04_Figures/F03/a_script/ring_helpers.R:43` assigns `padj = dep[["pi_score_<contrast>"]]`
  — eleven lines from a genuine fgsea BH q-value.

Everything downstream that calls a π-selected protein "significant" inherits this: F02's
entire shared-global / protein-specific-divergence thesis, F04's 55 divergent, F05's 86.

**Suggested fix.** Report BH as primary and state plainly that nothing survives. Keep
π < 0.05 as an explicitly-labelled exploratory ranking for feeding enrichment. Rename the
`padj` column. Add `limma::treat(fit, lfc = log2(1.2))` as the error-rate-controlled,
effect-size-floored analysis — it will likely return zero, which is the correct and far
stronger statement.

### B. F06 Panel B: fgsea with WGCNA modules as the gene sets

`panels/panel_B_nes.R` runs `fgsea::fgseaMultilevel()` using **the 12 WGCNA modules as
gene sets**, ranked by moderated-t from the same 45-sample matrix that defined those
modules. Committed q-values reach **1.9 × 10⁻⁴⁴**. Nowhere else in the pipeline does any q
drop below 0.24.

fgsea's competitive null assumes within-set gene independence. A WGCNA module is *by
construction* the maximal block of co-expressed proteins in that exact matrix — the
strongest possible within-set correlation. The p-value divides by a variance that is far
too small.

The **NES directions are probably genuine** (9 of 12 modules reverse between acute and
training — a real and interesting pattern). The *significance* is manufactured.

**Suggested fix.** `limma::fry` with the subject block and the `duplicateCorrelation`
consensus — already used correctly in `04_Figures/functions/concordance.R`.

### C. The reported prediction q is not the q

`prediction_shared/_leaf.R:104` applies `p.adjust(perm_p_q2, "BH")` **within one 36-cell
phase leaf**. Each phase is a separate `run.R`; nothing pools them.

| BH family | min q |
|---|---|
| One 36-cell phase leaf (**as reported**) | **0.358** |
| All 108 continuous cells | **0.806** |

The headline permutation p is **0.00995 = 2/201** — two permutations off the resolution
floor of a 200-permutation null. There is no max-statistic null anywhere, despite
`perm_matrix()` already reusing one seeded permutation set across cells, so the null q² of
every cell under each shuffle already exists. Building it is one `apply()`.

Also: only **7 of 108** continuous cells reach a positive Q² at all.

### D. The F06 eigengenes are transductive everywhere they are used

WGCNA is fit once on all 45 samples, then its eigengenes are used as features inside
leave-one-subject-out loops in F06 and in both 05_test prediction arms. The held-out
subject helped define its own features. No label leaks, so it is transduction rather than
leakage — but the Q²/AUC are optimistic and are **not out-of-sample in the strict sense**.
The eigengene cells are precisely the ones that top the tables.

**Suggested fix.** Refit WGCNA inside each fold, or drop the out-of-sample language.

### E. There is no contaminant filtering

No cRAP / `contaminants.fasta` removal and no decoy/reverse-hit step exists anywhere in
stage 01. Survivors in the final 1,896-protein matrix:

- **KRT1, KRT2, KRT10, KRT8** — KRT1/2/10 are the three canonical skin-handling
  contaminants at the top of every cRAP list. They survive because the blood rule only
  asks "secreted to blood / immunoglobulin / erythrocyte-high", and keratins are none of
  those. Keratin abundance in a biopsy is driven by sample handling, not biology.
- **HBG2** (fetal haemoglobin γ-2, HPA blood concentration 7 g/L). It slips the gate
  because the erythrocyte test keys on *adult*-erythrocyte scRNA (HBG2 nCPM = 0 — the
  fetal chain isn't transcribed in adult RBCs) rather than on the protein's own plasma
  abundance. HBB and HBA1 were correctly caught; HBG2 was not.

These columns enter the cycloess reference, all nine DEP contrasts, the WGCNA network as
legitimate nodes, and the ORA/GSEA universe.

### F. `HPA_annotations.tsv` has no provenance

`01_run_filtering.R:97-106` drops any protein absent from the HPA export — a 12,893-row
**partial** file with no recorded version, query string, or retrieval date, and no
`PROVENANCE.md`. It is load-bearing: 127 proteins dropped.

### G. Two degenerate cells in the class-prediction arm

Baseline protein and eigengene `glmnet` cells return **AUC exactly 0.000, CI 0.000–0.000,
perm p = 1.0**. A zero-width interval at exactly zero is a degenerate fit, not a
perfectly-inverted classifier. Worth a look at the harness.

### H. F02's composite tags do not match its panel filenames

`HRvLR_F02_composite.R` adds `pE` before `pD`, so patchwork's auto-tag **D** displays the
file named `panel_e_dep_counts.png` and tag **E** displays `panel_d_effectsize.png`.
Cosmetic, but confusing to anyone matching figure to file. The narrative documents the
mapping rather than renaming, to avoid touching figure code.

---

## Nine of sixteen rendered HTML files are stale

The `.qmd` sources are newer than their `.html` renders. Re-render before circulating.

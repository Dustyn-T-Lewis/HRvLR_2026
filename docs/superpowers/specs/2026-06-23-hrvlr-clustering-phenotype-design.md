# HRvLR clustering → phenotype linkage — design

Date: 2026-06-23
Stage: `05_Clustering/`
Status: design approved (roles + structure); two inputs pending before `d_integration`.
Constraint: pure R, no Python (MEFISTO/MOFA2 dropped — Python-only backend).

## Purpose

Group the HRvLR proteome into protein sets, summarise each set by an eigengene,
and relate those eigengenes to subject phenotype with a mixed model. The phenotype
linkage is the goal; the grouping is the means. Three grouping engines are run, but
only one carries the inferential claim.

## Design constraints (carry these everywhere)

- **Effective n = 16 subjects** (HR 8 / LR 8), not 45 samples. Two HR subjects are
  partial after QC (S28: T1 only; S29: T2/T3 only), so the 45 rows are mostly- but
  not-fully-paired repeated measures. Every method choice below is made for ~16
  subjects, not 45.
- **HR-only attrition.** All three dropped samples are HR (subjects S28, S29, high
  missingness). Dropout is QC-legitimate but not group-blind, so HR-specific results
  lean on the existing 48-vs-45 outlier-sensitivity fit
  (`03_DEP/.../02_dep_reports.R` → `11_outlier_sensitivity.csv`).
- **The whole stage is hypothesis-generating, not confirmatory.** At n = 15 this is
  the honest label for the figure and the manuscript text.
- **Pure R.** No reticulate/Python dependencies anywhere in the stage.

## Input and background

- **Clustering input set** = proteins with pi-score < 0.05 in *any* group (the
  pi-gated set). Pattern-agnostic, so it does not pre-bake the convergence pattern.
- **ORA background** = that *same* pi-gated set, for every engine. Not the 1,885
  quantified proteome, not the genome. One background, applied once, in
  `d_integration`.

## Architecture

Each engine is independent and phenotype-blind (except the quarantined supervised
one). Each writes the **same two artifacts** to its `c_data/`, which is the only
interface the downstream step depends on:

1. `membership.csv` — `protein_id, group_id, membership_weight` (weight = soft
   membership for Mfuzz, 1 for hard assignment).
2. `eigengene.csv` — `sample_id, subject, group_arm, timepoint, group_id, ME` —
   the per-sample eigengene (first PC of the group, sign-fixed; see below).

`d_integration/` consumes those, and nothing else, from all engines. This is what
keeps the LMM, the multiple-testing correction, the null, the concordance check, and
the ORA in one place under one rule instead of three diverging copies.

## Engines and roles

| Dir | Engine | Role | May claim |
|---|---|---|---|
| `a_wgcna_paired/` | WGCNA, paired design, within-subject centered | **PRIMARY — inferential** | "modules track phenotype"; BH across modules; permutation null; β + CI |
| `b_mfuzz_gap/` | Mfuzz on HR−LR gap | companion — temporal | convergence/divergence story; concordance only |
| `c_supervised/` | sPLS / DIABLO (mixOmics) | **QUARANTINE — exploratory** | hypothesis-generating only; nested CV; never the headline |

Why WGCNA as primary (after dropping MEFISTO): it is the canonical module→eigengene→
trait method, pure R, already installed. It does not model time inside the grouping —
that is handled downstream by within-subject centering (module definition) plus the
`timepoint` term and `(1|subject)` in the LMM. The honest cost: at n ≈ 15 WGCNA sits
at its stability floor, which is why the stage is labeled hypothesis-generating and
why every association is checked against the permutation null. Anchor for the paired
handling: Li et al. 2018 (*Sci Rep*).

## Per-engine defaults

**WGCNA paired (`a_wgcna_paired/`) — PRIMARY**
- **Within-subject centering first** (subtract each subject's across-timepoint mean
  from the abundance matrix). This is the step that stops modules from encoding
  between-subject identity — the dominant failure mode at this n.
- Signed network; soft-threshold power at scale-free topology R² > 0.8;
  `minModuleSize` 20; `mergeCutHeight` 0.25; `deepSplit` 2.
- Eigengene = `moduleEigengenes()$eigengenes`, sign-fixed so positive ME = higher mean
  abundance of the module's positive-loading proteins.
- Report module-count sensitivity to power and `mergeCutHeight`.
- Alternative / switch: if modules are unstable across power choices, raise
  `minModuleSize` and report fewer, larger modules rather than chasing many small ones.

**Mfuzz on HR−LR gap (`b_mfuzz_gap/`) — companion**
- Feature per protein = 3-point HR−LR contrast trajectory `[ΔT1, ΔT2, ΔT3]` built
  from the limma moderated logFCs already computed in `03_DEP` (denoised,
  design-aware, sidesteps the missForest-manufactured-structure risk).
- z-score per protein; `c` bounded to ~4–6 by the 3-point geometry, cross-checked
  with `Dmin`; `m` via `mestimate`; membership cutoff ≥ 0.5 (reported); seed-stability
  across ≥ 25 seeds.
- Also produce the group-mean 6-point profile clustering as the *descriptive* shape
  figure (illustration only — the gap trajectory makes the convergence claim).
- Reporting rule (pre-committed): gap = inference-eligible within this companion;
  shape = illustration. Do not swap based on which looks nicer.

**Supervised (`c_supervised/`) — quarantine**
- sPLS / DIABLO with nested cross-validation; report CV error honestly.
- Labeled exploratory throughout. Used to generate candidate proteins, not to assert
  phenotype association. At n = 15 these overfit; nested CV is the only defense and is
  itself unstable here.

## `d_integration/` — the shared, honest downstream

- **LMM** per group (primary engine first): `ME ~ phenotype + group + timepoint +
  (1|subject)`, `lmerTest` for inference, standardised β + partial/semi-partial R²
  for effect size. Random intercept only — n = 15 will not support random slopes.
- **Multiple testing**: BH across the primary engine's modules. Companion engines are
  *not* added to the same correction — they are reported as concordant/discordant with
  the primary, not as extra tests.
- **Permutation null (falsification test)**: permute phenotype across the 16 subjects
  (subject-level, preserving each subject's actual timepoint block — two HR subjects
  have fewer than 3), recompute the largest
  |module–phenotype| association, B = 1000. The observed top association must beat this
  null or the "tracks phenotype" story is not supported.
- **Cross-engine concordance**: cross-tab membership across engines (adjusted Rand
  index, Jaccard on group overlaps). If WGCNA modules and Mfuzz clusters are
  near-identical, say so and collapse to one figure — do not present duplicate panels.
- **ORA**: `clusterProfiler::enricher` against the pi-gated background, GO BP/CC/MF +
  Reactome, BH within each gene-set source, **no `clusterProfiler::simplify()`**.

## Open inputs to confirm before `d_integration` runs

1. **Phenotype variables**: which trait(s), and measured per-timepoint (T1/T2/T3) or
   once per subject. This sets whether phenotype is a within-subject term or a
   subject-level constant, and it sets the permutation scheme.
   *Recommended default*: one primary trait, subject-level if that is how it was
   measured; add a small pre-named panel only with BH across traits.
2. **Timepoint coding** for the LMM `timepoint` term and the Mfuzz gap: T1 baseline,
   T2 post-training, T3 acute. T3 is acute response, not simply "later than T2."
   *Recommended default*: treat T1 → T2 as the training axis and T3 as a separate
   acute level, not an ordered third time. (This mattered most for MEFISTO's smooth
   covariate, now dropped; it still affects how `timepoint` is coded downstream.)

## Decisions ranked by conclusion-impact

1. Effective n = 16, unbalanced (sets everything; forces hypothesis-generating framing).
2. Phenotype binding + permutation scheme (open input #1).
3. Within-subject centering for WGCNA (decides whether the *primary* modules are
   signal or subject identity).
4. Which engine is inferential vs companion vs quarantine (locked: WGCNA primary).
5. Trajectory representation for Mfuzz = HR−LR gap from limma logFCs (decides what
   "convergence" means and dodges the imputation artifact).
6. ORA background = pi-gated set (decides enrichment validity).
7. Timepoint coding / T3 handling (open input #2).
8. Per-engine tuning (power, `mergeCutHeight`, `c`, `m`) — real but lower-impact.

## Fragilities and anchors

- **Pseudo-replication** (3 timepoints/subject): one network across 45 samples
  inflates apparent module stability; handled by within-subject centering (WGCNA) and
  by `(1|subject)` in the LMM. Anchors: Langfelder & Horvath 2008; Li et al. 2018.
- **Circularity / double-dipping**: phenotype-blind grouping → eigengene → LMM is
  *not* double-dipping (selection is phenotype-independent; Kriegeskorte et al. 2009
  does not apply). Supervised covariation *is*, which is why it is quarantined.
- **Imputation-manufactured structure**: missForest borrows across proteins, so
  clustering raw imputed values can fabricate co-clustering. Mfuzz uses limma
  estimates instead of raw imputed values; high-missingness proteins are flagged.
- **3-timepoint geometry**: fuzzy c-means was built for long series; on 3 points the
  shape space is nearly enumerable, so `c` is geometry-bounded. Mfuzz is kept for the
  soft-membership down-weighting, which matters at this n. Anchor: Kumar & Futschik 2007.
- **WGCNA at the n-floor**: ~16 independent units is right at the documented minimum;
  modules are reported with stability checks and never claimed beyond the permutation
  null. (Empirically the 5 modules were invariant across powers and mergeCutHeights.)

## What would falsify the convergence/divergence story

The permutation null above. If the observed top module–phenotype association does not
exceed the subject-permuted null, and if the HR−LR gap clusters do not separate from
a label-permuted gap, the convergence/divergence reading is not supported and is
reported as null.

## Deliverable / figure sequence (only after the design holds)

1. Input + roster panel: pi-gated set size, n = 15, HR-only attrition caveat.
2. Primary (WGCNA): modules, soft-threshold/scale-free fit, module–phenotype
   associations with the permutation null, hub proteins.
3. Companion temporal (Mfuzz gap): convergence/divergence clusters + descriptive
   shape profiles.
4. Concordance: cross-engine membership cross-tab (justifies how many engines survive
   into the final figure).
5. ORA on the surviving groups (consistent background).
6. Supplement: supervised exploratory result, labeled non-confirmatory.

## Directory layout

```
05_Clustering/
├── a_wgcna_paired/   {a_script, b_reports, c_data}   # PRIMARY
├── b_mfuzz_gap/      {a_script, b_reports, c_data}   # companion: temporal shape
├── c_supervised/     {a_script, b_reports, c_data}   # QUARANTINE: exploratory
└── d_integration/    {a_script, b_reports, c_data}   # shared: LMM, null, BH,
                                                       #   concordance, ORA
```

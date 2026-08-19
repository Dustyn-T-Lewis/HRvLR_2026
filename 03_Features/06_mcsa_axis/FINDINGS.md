# mCSA axis — findings

Run 2026-08-19 under the rules fixed in PREREG.md. One figure carries
everything: `b_reports/F_mcsa_axis.{png,pdf}` with its legend beside it.
Seed 42, 931 complete-case proteins, blood index in every model.

## The phenotype is two axes, and one subject sets how far apart

Whole-muscle CSA change and the fibre pair are nearly independent over the
16 subjects (Spearman 0.09 with Δ fCSA I, 0.20 with Δ fCSA II), and the two
strength tests share no axis with anything or each other (−0.39).

That near-independence rests heavily on LR_S14. Removing it moves Δ mCSA
against Δ fCSA I from 0.09 to 0.32, against Δ fCSA II from 0.20 to 0.45, and
against the composite from 0.49 to 0.69 — the three largest shifts of any
pair in the matrix (`c_data/01_phenotype_geometry.csv`). LR_S14 posts the
study's largest whole-muscle gain (z = +1.63) with near-worst fibre change
(z = −1.75, −1.80). Fibres down and whole muscle up in one leg reads as a
measurement discrepancy at least as readily as two biological axes.

## Δ mCSA is an ingredient of the arm definition, not a phenotype outside it

Regressing `comp_hypertrophy` on the z-scored CSA measures returns R² = 0.964
with weights 4.98 (fCSA I), 2.41 (fCSA II) and 4.12 (mCSA); dropping the
whole-muscle term takes R² to 0.737. Δ mCSA alone separates HR from LR at
AUC 0.83, Wilcoxon p = 0.028, its top eight holding six of the eight real HRs
(`c_data/01_arm_separation.csv`).

So the reconciliation that motivated this stage — that F04's null and F06's
Δ mCSA leads concern different phenotypes — does not hold as stated. What
holds is narrower: Δ mCSA carries 36% of the composite's formula and 12% of
the fitted factor's weighting (`c_data/01_composite_weights.csv`). It is a
third of the recipe and almost none of the ordering.

## Proteins against Δ mCSA: one survivor, concurrent, at T2

Eight cells, four configs by two subject sets, BH within a config across the
931 (`c_data/02_scan_summary.csv`):

| Config | n | min BH q | nominal p < 0.05 | survivors |
| --- | --- | --- | --- | --- |
| T1 | 15 | 0.944 | 31 | 0 |
| T2 | 15 | **0.045** | 103 | **1** |
| T3 | 15 | 0.246 | 112 | 0 |
| Δ (T2 − T1) | 14 | 0.992 | 25 | 0 |
| T1, LR_S14 dropped | 14 | 0.999 | 16 | 0 |
| T2, LR_S14 dropped | 14 | 0.506 | 60 | 0 |
| T3, LR_S14 dropped | 14 | 0.420 | 99 | 0 |
| Δ, LR_S14 dropped | 13 | 0.936 | 36 | 0 |

Chance predicts 47 nominal hits per cell. T1 and the delta return fewer than
chance in both subject sets; T2 and T3 return more.

The survivor is **PSME1** (Q06323, proteasome activator complex subunit 1) at
T2: slope 0.225 log2 per unit Δ mCSA, t = 5.48, p = 4.8 × 10⁻⁵, BH q = 0.045.
It cleared the second gate as well — 0 of 200 subject-permutations reached its
|t|, empirical p = 0.005 — so the permutation armed and it passed.

What it survives (`c_data/02_survivor_checks.csv`):

- Leave-one-subject-out: t between 4.47 and 8.86 across all 15 removals. No
  single subject carries it.
- Within each arm: Spearman 0.86 in HR (n = 7), 0.98 in LR (n = 8).
- Adding the arm label: the arm term is not significant (p = 0.27) and the
  Δ mCSA slope rises to 0.248 (p = 7.7 × 10⁻⁵). This is not the HR-vs-LR
  contrast in another coat, which matters because Δ mCSA separates the arms
  at AUC 0.83.

What it does not survive, and what it never claimed:

- **It is not a forecast.** At T1 PSME1 ranks 242nd of 931 (p = 0.25,
  q = 0.944). The same protein and the same subjects carry no baseline signal.
  Δ mCSA is measured over T1 to T2, so only T1 precedes it.
- **It does not clear BH with LR_S14 removed** (q = 0.506), though it stays
  the top-ranked protein in that cell (p = 5.4 × 10⁻⁴, slope 0.196 against
  0.225). The effect size barely moves; the threshold does.
- **It is one survivor across eight cells.** BH ran within a config as
  pre-registered, not across the eight. Read against the eight-cell family it
  is a single 0.045.
- It does not replicate at T3 (rank 6, q = 0.372) or on the training delta
  (rank 33, q = 0.992), although the T3 sign and rank agree.

## What this changes, and what it does not

The fibre axis still has no proteome correlate anywhere tested: 2 leads in 396
F06 fibre-axis cells (0.5%, under the 5% chance rate), zero BH survivors
against the galamm latent factor (minimum q = 0.117), and zero at protein and
module level in the nine arm × timepoint contrasts. Splitting the phenotype
into two outcomes changed no F04 or F06 conclusion, because F04 fits no
continuous phenotype and F06 never pooled the six traits.

What the split did buy is one protein. PSME1's abundance 72 h after training
tracks how much whole-muscle area a subject gained, within arm, with the blood
index partialled out, at a level that clears BH and its own permutation. It is
a concurrent correlate on a phenotype whose own measurement is in question at
one subject, drawn from eight cells. It is worth a sentence and a supplement,
not a claim about responders.

## Deviations from PREREG.md

- `c_data/02_survivor_checks.csv` and `survivor_checks()` were added after a
  survivor appeared. They were not pre-declared. They diagnose a hit rather
  than gate one, and no gate moved because of them.
- The figure's third panel was pre-declared as an observed-versus-permutation
  envelope. Permutation ran only behind the BH gate, so no envelope exists
  for the 931; the panel shows observed against uniform expected quantiles
  instead, which needs no extra fitting and is labelled as such.
- The fourth panel was pre-declared as the F06 lead rate by outcome. With a
  survivor to show, that panel became PSME1 across the three timepoints and
  the F06 rates moved into the legend. Panel count and the four-panel limit
  are unchanged.
- Nothing else: subset, covariate, configs, thresholds, gates and stop rules
  ran as written.

## Proposed README paragraph (not inserted)

> Splitting the phenotype into a fibre axis and whole-muscle CSA changes no
> conclusion in F04, which fits no continuous phenotype, or in F06, which
> never pooled the six traits (`03_Features/06_mcsa_axis/`, pre-registered).
> Whole-muscle CSA change is an ingredient of the composite the arms were
> split on rather than a phenotype outside it: it carries 36% of the
> composite's formula, separates the arms at AUC 0.83, and contributes 12% of
> the fitted factor's weighting. Scanning the 931 complete-case proteins
> against it with the blood index partialled out returns one BH survivor in
> eight cells — PSME1 at T2, q = 0.045, permutation p = 0.005, holding under
> leave-one-out and within both arms, and not explained by the arm label. At
> T1 the same protein ranks 242nd (p = 0.25), so it is a concurrent correlate
> of whole-muscle growth, not a baseline forecast of it. The fibre axis has no
> protein correlate at any level tested.

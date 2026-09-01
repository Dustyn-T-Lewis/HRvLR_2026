**The mCSA axis: what the phenotype is made of, and what the proteome tracks.**

**A.** Each subject's fibre axis (mean z of Δ fCSA I and Δ fCSA II) against
whole-muscle CSA change, coloured by arm. Spearman ρ = 0.12 over all 16
subjects and 0.36 with LR_S14 removed. LR_S14 posts the study's largest
whole-muscle gain (z = +1.63) alongside near-worst fibre change (z = −1.75,
−1.80); fibres down and whole muscle up in the same leg is a measurement
discrepancy as readily as two biological axes, and the three correlations that
move most when any one subject leaves are all Δ mCSA pairs.

**B.** The three CSA measures as shares of two weightings. Left, the
composite's own formula, recovered by regressing `comp_hypertrophy` on the
z-scored measures (R² = 0.964; 0.737 without the whole-muscle term). Right,
the galamm measurement model's loadings, normalised across the same three
(`03_Features/04_galamm_pilot`). Whole-muscle CSA carries 36% of the formula
and 12% of the fitted factor. It is an ingredient of the composite the HR/LR
arms were split on, not a phenotype outside it: on its own it separates the
arms at AUC 0.83, Wilcoxon p = 0.028.

**C.** Quantile plot of the 931 complete-case proteins regressed on `d_mcsa`
with the blood index partialled out, one curve per config, dashed with LR_S14
removed. Expected quantiles are uniform, not a permutation envelope. Only the
T2 curve crosses BH q = 0.05. Δ `d_mcsa` spans T1 to T2, so T1 is the only
config that precedes its outcome and the only one that could forecast it; its
curve sits on the diagonal.

**D.** PSME1 against `d_mcsa` at each timepoint, same subjects throughout, BH q
within that config across the 931. Grey band is the ordinary least-squares
interval.

**What the survivor is and is not.** PSME1 clears BH at T2 (q = 0.045) and
beats its subject-permutation null (0 of 200 permutations reached |t| = 5.48,
empirical p = 0.005). It holds when any one subject is removed (leave-one-out
t 4.47 to 8.86), inside both arms separately (ρ = 0.86 HR, 0.98 LR), and when
the arm label is added to the model, where the arm term is not significant
(p = 0.27) and the `d_mcsa` slope rises from 0.225 to 0.248. It is therefore
not the HR-vs-LR contrast in another form.

It is one survivor across eight cells, BH applied within a config as
pre-registered rather than across the eight. It does not clear BH with LR_S14
removed (q = 0.506) although it stays the top-ranked protein there
(p = 5.4 × 10⁻⁴, slope 0.196). At T1 it ranks 242nd of 931 (p = 0.25), so it
is a concurrent correlate of whole-muscle growth, not a baseline forecast of
it.

**Context.** The fibre axis has no proteome correlate at any level tested: the
F06 prediction screen returns 2 leads in 396 fibre-axis cells (0.5%, below the
5% chance rate) against 32 of 132 for `d_mcsa`; the galamm pilot's latent
factor, which is the fibre axis, associates with none of the 931 proteins
(minimum BH q = 0.117); and the nine arm × timepoint contrasts have zero BH
survivors at protein and module level.

Scripts: `03_Features/06_mcsa_axis/a_script/`. Pre-registration: `PREREG.md`.
Data: `c_data/01_*.csv`, `c_data/02_*.csv`.

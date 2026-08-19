# 05_phenotype_modules — findings

Run 2026-08-19. Label-free follow-up to the galamm pilot: can protein
clusters or modules, built with no HR/LR label, associate with the
continuous phenotype outcomes, tested independently at each timepoint?
Figure: `b_reports/F_phenotype_modules.{png,pdf}`. 14 of 16 subjects
carry all three timepoints; B = 1000 subject-permutation nulls; seed 42.
The consistency call (same sign at T1, T2 and T3, empirical p < 0.05 at
two or more) runs identically inside every permuted cohort, so the
count of consistent features has its own null. Timepoints share
subjects: consistency is stability of an association, not independent
replication.

## Module level: nothing beats the consistency null

One module is consistent for one phenotype — brown vs d_mcsa (rho 0.57,
0.74, 0.67 at T1, T2, T3) — against a null median of 0, consistency
p = 0.057. Two caveats keep it a near-miss rather than a result: it
does not clear 0.05, and brown is one of the ten modules that fail
leave-one-subject-out network rebuilding (only pink and turquoise
survive, and neither associates with anything). No other module is
consistent for any phenotype.

## Pathway level: the counts sit inside their nulls

Of 57 Hallmark singscore sets x 6 phenotypes, five sets reach p < 0.05
twice for one phenotype: Interferon Gamma Response and E2F Targets
(positive, d_mcsa, consistency p = 0.09), Peroxisome (negative,
composite, p = 0.27), Bile Acid Metabolism and Adipogenesis (negative,
d_1rm_ext, p = 0.10). Every count is inside its permuted distribution.
The d_mcsa concentration echoes F06's — and like F06's survivors these
are concurrent correlates, not baseline forecasts: the T1-only column
carries none of the boxed cells on its own.

## Delta clusters: a clean negative

Clustering the 931 complete-case T2 - T1 deltas splits the proteome
into an up-mover and a down-mover cluster (k = 2 by silhouette,
n = 557 and 374). Neither cluster's per-subject mean delta reaches
permutation p < 0.05 for any phenotype. Proteins that move together
over training do not move with how much anyone grew.

## Bottom line

Removing the labels does not rescue the proteome-phenotype link. The
only signal shape that recurs is the one F06 already found — d_mcsa,
concurrent, modest — and here it peaks at consistency p = 0.057 on a
module that does not survive resampling. The null stands at a fourth
level of description (proteins, latent factor, modules/pathways, delta
clusters).

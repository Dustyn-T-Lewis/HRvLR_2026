# Proteome

If people are grouped by their proteome alone, what else differs between the groups?

## What was clustered

The 15 module eigengenes from `01_Processing/06_Module_Building`, with no phenotype or group
label, five times over:

| Input | Each person's values |
|---|---|
| `T1`, `T2`, `T3` | eigengenes at that time point |
| `training` | eigengenes at T2 minus T1 |
| `acute` | eigengenes at T3 minus T2 |

The gap statistic selects K = 1 at every input. K = 2 is imposed at every input with the recorded
reason: so the exploratory folder can compare the two halves the proteome suggests; read each split
with its bootstrap Jaccard. The two clusters are `C1` and `C2`, numbered by PAM, so the names
carry no direction and do not match across inputs.

The method (PAM, gap statistic, bootstrap Jaccard) is in [`../README.md`](../README.md) and
[`01_Clustering`](01_Clustering/README.md).

## Validity: circular, caution, clean

Here the groups come from the proteome, and the later stages test the proteome. Whether a
comparison means anything depends on which time points it shares with the clustering.
`validity()` in [`R/inputs.R`](R/README.md) tags every comparison against the input it was run on:

| Tag | Rule | Why |
|---|---|---|
| circular | the comparison uses exactly the clustering's time points | The clustering sorted people on those values, so the clusters differ there by construction, even in pure noise. The clustering made the difference; the result describes the split and tests nothing. |
| caution | it shares a time point with the clustering but is not the same comparison | Regression to the mean. The clustering sorted people partly on chance high or low values at the shared time point. That chance enters the comparison too, and at an unshared time point it does not repeat, so the clusters can seem to differ, or to change differently, with no real difference. |
| clean | it shares none | The noise the clustering sorted on does not enter the comparison. A difference here is a test of whether the split carries over to another time point or window. |

`uses()` reads the time points from the row's `contrast`, `task` or `window` name:

| Uses | Names |
|---|---|
| T1 | `Baseline_C1vC2` (contrast and classification task); the `baseline` and `baseline_level` windows; PERMANOVA window `T1` |
| T2 | `Trained_C1vC2` (contrast and task); `post_level`; PERMANOVA window `T2` |
| T3 | `Acute_C1vC2` (contrast and task); `post_bout_level`; PERMANOVA window `T3` |
| training (T1, T2) | `Training_C1`, `Training_C2`, `Training_Interaction`; the `Training_change_C1vC2` task; the `training` windows |
| acute (T2, T3) | `Acute_C1`, `Acute_C2`, `Acute_Interaction`; the `Acute_change_C1vC2` task; the `acute` windows |
| all (T1, T2, T3) | PERMANOVA windows `all` and `subject mean` |

Against each input that gives:

| Clustered at | T1 | T2 | T3 | training | acute | all |
|---|---|---|---|---|---|---|
| `T1` | circular | clean | clean | caution | clean | caution |
| `T2` | clean | circular | clean | caution | caution | caution |
| `T3` | clean | clean | circular | clean | caution | caution |
| `training` | caution | caution | clean | circular | caution | caution |
| `acute` | clean | caution | caution | caution | circular | caution |

Read circular rows as description, caution rows with the regression-to-the-mean risk in mind, and
clean rows as the tests. All of them inherit the imposed K: the clusters are halves of one spread,
to be read with their Jaccard, so even a clean result is exploratory.

## Per-input folders and stacked workbooks

Steps `02_Design_Contrasts` to `07_Trait_Associations` call the shared function once per input,
each into its own folder (`03_Protein_Contrasts/T1/c_data/protein_contrasts.xlsx`). `stack_inputs()`
in [`R/inputs.R`](R/README.md) then joins the five workbooks into
`<step>/c_data/<name>_all_inputs.xlsx`: every result sheet with `input` first and, where a row
names a contrast, task or window, its `validity`. `07_Trait_Associations` leaves `within_group`
out of the stack, because five inputs pass xlsx's row limit; each input's workbook keeps it.

## Steps

| Step | Asks |
|---|---|
| [`01_Clustering`](01_Clustering/README.md) | how people split at each input |
| [`02_Design_Contrasts`](02_Design_Contrasts/README.md) | which comparisons |
| [`03_Protein_Contrasts`](03_Protein_Contrasts/README.md) | which proteins differ |
| [`04_Set_Contrasts`](04_Set_Contrasts/README.md) | which gene sets shift |
| [`05_Module_Contrasts`](05_Module_Contrasts/README.md) | which modules shift |
| [`06_Classification`](06_Classification/README.md) | can one protein, set or module tell C1 from C2 |
| [`07_Trait_Associations`](07_Trait_Associations/README.md) | does a feature's link to a trait differ by cluster |
| [`08_Cluster_Phenotype`](08_Cluster_Phenotype/README.md) | do the clusters differ in size or strength |
| [`Figures/F1_Overview`](Figures/F1_Overview/README.md) | the folder on a page, once per input |

The outputs of this folder are not tracked in git (`.gitignore`); `Rscript run_all.R 03` rebuilds
them with the two muscle folders.

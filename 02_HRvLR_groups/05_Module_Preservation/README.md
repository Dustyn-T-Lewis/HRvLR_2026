# 05_Module_Preservation

Does each group's co-expression structure reappear in the other group?

## The idea

01_Processing builds one network from the whole cohort. This stage builds a network inside each group
separately, with the same settings, and asks whether a module found in one group is still a
module in the other. A module that is missing in one group would point to proteins that are
coordinated in one kind of responder and not the other.

A co-expression network starts from the correlation between every pair of proteins across
biopsies. Here the correlation is the biweight midcorrelation (bicor), which down-weights outlying
points, on subject-centred abundance (each person's own mean removed, so the network reflects how
proteins move together within people). People with only one biopsy centre to zeros and are left
out. The correlations are raised to a soft-thresholding power so strong links dominate and weak
ones fade; the power is chosen by `WGCNA::pickSoftThreshold()`, with the WGCNA FAQ's table for
signed networks as the fallback when no power fits. Proteins are then clustered into modules and
named by colour.

`WGCNA::modulePreservation()` takes one group's modules (the reference) and looks at the same
proteins in the other group's data (the test). It asks two things: are the module's proteins still
more correlated with each other than with the rest (density), and do they keep the same pattern of
connections (connectivity)? Each statistic is set against 200 permutations of the module labels,
which give the values a random set of proteins of the same size would reach. A Z score is the
observed value minus the permutation mean, divided by the permutation standard deviation.
Zsummary combines the density and connectivity Z scores. The bands follow Langfelder et al.
(2011):

| Zsummary | Reading |
|---|---|
| above 10 | strong preservation |
| 2 to 10 | weak to moderate |
| below 2 | no evidence of preservation |

Zsummary grows with module size, so a large module scores higher than a small one equally well
preserved. The median rank ranks the modules against each other on the observed statistics and
does not depend on size; lower means better preserved.

Each group network rests on 20 to 24 samples. That is little power to show a module is absent.
A Zsummary above 2 says a module is found again in the other group; a Zsummary below 2 is weak
evidence that it is missing; neither shows the two networks are the same.

Module colours are assigned afresh in each group's network, so HR's blue is not LR's blue and
neither is the 01_Processing blue. Each group module is labelled with the full-cohort module it overlaps
most, so the names can be read together.

## What the code does

| Step | Function | Detail |
|---|---|---|
| Group data | `apply()`, `ave()` | Subject-centred abundance per group, from the blood-corrected complete matrix; people with one biopsy left out. |
| Soft power | `WGCNA::pickSoftThreshold(powerVector = 1:20, networkType = "signed", corFnc = "bicor")` | Falls back to 18, 16, 14 or 12 by sample count when no power reaches the fit criterion. |
| Modules | `WGCNA::blockwiseModules()` | Signed network and signed TOM, bicor (`maxPOutliers = 0.05`), `deepSplit = 2`, `minModuleSize = 30`, `mergeCutHeight = 0.15`, PAM stage on, one block, `randomSeed = 54321`: the 01_Processing settings. |
| Preservation | `WGCNA::modulePreservation(referenceNetworks = 1:2, nPermutations = 200, randomSeed = 42)` | Both directions: HR modules tested in LR, LR modules tested in HR. grey (unassigned) and gold (WGCNA's random reference module) are dropped. |
| Bands | `cut(z_summary, c(-Inf, 2, 10, Inf))` | none, weak to moderate, strong. |
| Naming | `dplyr::count()`, `slice_max()` | Each group module's best-overlapping full-cohort module and the share of its proteins there. |

## How to read the outputs

- Read Zsummary against module size. A small module below 2 may be missing or may simply be too
  small to show at this sample size; the median rank helps separate the two.
- A module strongly preserved in both directions is shared structure. A module preserved in one
  direction only says that group's module is not recovered as a unit in the other, which can mean
  it is split across several of the other group's modules.
- The test compares correlation structure, not abundance. A module can be fully preserved while
  its level differs between the groups; that is the question `04_Module_Contrasts` asks.
- This stage feeds no figure. It is a check on whether one network can describe both groups.

## Files

| | |
|---|---|
| Reads | `01_Processing/03_Quantification/c_data/DAList_corrected.rds` (blood-corrected complete matrix), `01_Processing/06_Module_Building/c_data/modules.rds` (sample metadata, full-cohort membership) |
| Writes | `c_data/module_preservation.xlsx`, `b_reports/module_preservation_report.pdf` |
| Run | `Rscript 02_HRvLR_groups/05_Module_Preservation/a_script/module_preservation.R` |

```
05_Module_Preservation/
├── a_script/module_preservation.R
├── b_reports/module_preservation_report.pdf
└── c_data/module_preservation.xlsx
```

## Workbook: `module_preservation.xlsx`

`read_me` lists the sheets; `columns` defines every column; `input_manifest` and
`package_versions` close the workbook.

| Sheet | One row is | Key columns |
|---|---|---|
| `preservation` | a group module tested in the other group | `reference`, `test`, `module`, `n_proteins`, `z_summary`, `median_rank`, `preservation` (band), `full_module`, `full_share` |
| `arm_networks` | a group's network | `arm`, `n_samples`, `power`, `n_modules` |
| `arm_membership` | a protein in a group's network | `arm`, `uniprot_id`, `module` |

## Report: `module_preservation_report.pdf`

1. One page: Zsummary (top) and median rank (bottom) against module size on a log scale, one
   column per direction. Each point is labelled with its module and, in brackets, the
   best-overlapping full-cohort module. Dashed lines mark Zsummary 2 and 10.

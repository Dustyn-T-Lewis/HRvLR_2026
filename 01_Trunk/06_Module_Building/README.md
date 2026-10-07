# 06_Module_Building

Which proteins rise and fall together within a person, and how does each such group move in each
biopsy?

## The idea

Gene sets come from databases. Modules come from the data: groups of proteins whose abundances
move together across the biopsies. WGCNA (weighted gene co-expression network analysis;
Langfelder & Horvath 2008) finds them in four steps:

1. Correlate every pair of proteins across samples. Here the correlation is biweight
   midcorrelation (bicor), a robust version of Pearson's that a single extreme sample cannot drive.
2. Turn correlations into connection strengths. The network is signed: a strong positive
   correlation is a strong connection, a strong negative one is none. Strengths are raised to a
   soft power, which shrinks weak correlations towards zero and keeps strong ones. The power is
   the lowest one at which the network looks scale-free (a few highly connected proteins, many
   sparsely connected), judged by a fit R² of 0.85.
3. Measure topological overlap (TOM): two proteins are close when they share many of the same
   neighbours, not only when they correlate with each other. This steadies the clusters.
4. Cluster proteins on TOM and cut the tree into modules. Each module gets a colour name; `grey`
   holds the proteins that fit no module. Modules whose summaries correlate closely are merged.

Two summaries describe a module:

- The eigengene is the first principal component of the module's proteins: one number per biopsy
  that best represents how the whole module moves. Its sign is aligned so it rises when the
  module's proteins on average rise.
- kME (module membership) is each protein's correlation with its module's eigengene. A protein
  with kME near 1 tracks the module closely; the highest-kME proteins are its hubs.

WGCNA treats samples as independent. These 45 biopsies come from 16 people, up to three each, and
subject identity drives the leading patterns in raw abundance. Fed raw abundance, WGCNA would
build modules that mostly say who a biopsy came from. So the modules are defined on
subject-centred abundance: each person's mean is subtracted from each of their biopsies, leaving
only how a protein moved within that person. A subject with one biopsy centres to zeros and is
left out of the definition.

The eigengenes are then scored on uncentred abundance across all 45 biopsies, so differences
between people, including between arms, stay testable in the arcs. kME is read on the
subject-centred data that defined the modules.

The input is limpa's complete matrix with the blood slope removed. Construction never sees a group
label or a phenotype.

## What the code does

| Step | Function | Detail |
|---|---|---|
| Centre | `ave()` | Each protein centred within subject, over subjects with more than one biopsy. |
| Soft power | `WGCNA::pickSoftThreshold()` | Powers 1 to 20, signed, bicor (`maxPOutliers = 0.05`). Uses the package's estimate, the lowest power whose scale-free R² clears 0.85; if none does, the WGCNA FAQ's signed power for the sample count (18 under 20 samples, 16 for 20-30, 14 for 31-40, 12 over 40). |
| Modules | `WGCNA::blockwiseModules()` | Signed network and signed TOM, bicor with `maxPOutliers = 0.05` and `pearsonFallback = "individual"`, `deepSplit = 2`, `minModuleSize = 30`, `mergeCutHeight = 0.15`, `pamStage = TRUE`, `maxBlockSize = 2500` (one block), `randomSeed = 54321`. Other settings are package defaults. |
| Eigengenes | `WGCNA::moduleEigengenes(excludeGrey = TRUE)` | On uncentred abundance, all 45 biopsies. |
| Membership | `WGCNA::moduleEigengenes()`, `WGCNA::signedKME()` | Eigengenes on the centred matrix, then each protein's bicor with its own module's eigengene. Grey proteins get no kME. |
| Subject dependence | `lme4::lmer(eigengene ~ timepoint + (1 \| subject))` | The intraclass correlation (ICC) of each scored eigengene: between-subject variance over between-subject plus residual variance, with time point fitted so its effect does not count as noise. |

## How to read the outputs

- Module colours are names, not ranks, and WGCNA reassigns them on every build. When the trunk
  changes, `turquoise` may hold different proteins; compare modules by membership, not by colour.
- A high ICC means an eigengene still mostly separates people, even though the module was defined
  within people. That is expected for a module scored on uncentred abundance and is the reason the
  arcs block on subject.
- kME is read on centred data and the eigengene on uncentred data. A protein's kME describes how it
  co-varies within a person; it is not its correlation with the eigengene values in `eigengenes`.
- `soft_threshold` shows the fit at every power; `chosen` marks the power used. The FAQ fallback
  applies only when no power's `SFT.R.sq` exceeds 0.85.
- What each module is about is answered in `07_Module_Characterisation`.

## Files

| | |
|---|---|
| Reads | `01_Trunk/03_Quantification/c_data/DAList_corrected.rds` |
| Writes | `c_data/modules.rds`, `c_data/module_building.xlsx`, `b_reports/module_building_report.pdf` |
| Run | `Rscript 01_Trunk/06_Module_Building/a_script/module_building.R` |

`modules.rds` holds `eigengenes` (module by sample matrix, grey excluded), `membership` (every
protein's module and kME), `meta` (sample sheet: `sample_id`, `subject`, `arm`, `timepoint`) and
`module_summary`. It is read by `07_Module_Characterisation`, by the module contrasts,
preservation, classification, trait associations and F2 figure in `02_HRvLR_groups`, by the module
steps and profile clusters in `04_Continuous_no_group`, and by each
`03_Cluster_groups/*/01_Clustering` (`Proteome` clusters people on the eigengenes).

```
06_Module_Building/
├── a_script/module_building.R
├── b_reports/module_building_report.pdf
└── c_data/
    ├── module_building.xlsx
    └── modules.rds
```

## Workbook: `module_building.xlsx`

`read_me` lists the sheets; `columns` defines every column; `input_manifest` and
`package_versions` close it.

| Sheet | One row is | Key columns |
|---|---|---|
| `module_summary` | a module (grey excluded) | `n_proteins`, `icc` |
| `membership` | a protein | `uniprot_id`, `gene`, `module`, `kme` |
| `eigengenes` | a sample | `sample_id`, then one eigengene per module |
| `soft_threshold` | a soft power | `SFT.R.sq`, `slope`, `signed_r2`, `mean.k.`, `chosen` |

## Report: `module_building_report.pdf`

1. Soft-threshold choice: signed scale-free R² and mean connectivity at each power, the 0.85 line
   and the chosen power.
2. Module sizes, each bar labelled with its eigengene's subject ICC.
3. Eigengene trajectories: one thin line per subject, arm mean and standard error at each time
   point, one panel per module.

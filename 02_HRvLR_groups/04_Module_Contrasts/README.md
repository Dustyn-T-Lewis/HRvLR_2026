# 04_Module_Contrasts

Which co-expression modules shift in each of the nine contrasts, and do their most central
proteins move most?

## The idea

A module is a group of proteins whose abundances rise and fall together across biopsies. The
modules come from WGCNA in `01_Trunk/06_Module_Building`, built without any group label or trait
and named by colour; grey holds the proteins that join no module and is left out here. Modules
were defined on subject-centred abundance (each person's own mean removed), so they capture
proteins that move together within people rather than proteins that differ between people.

The eigengene is a module's summary: the first principal component of its members, one number per
biopsy, signed to follow the members' average. It is scored on the blood-corrected abundance, not
the centred one, so differences between groups stay in it and can be tested. It is on its own
standardised scale, so an eigengene logFC is a difference in eigengene units, not a log2 fold
change.

kME (module membership) is the correlation between one protein and its module's eigengene. Near
1, the protein is central and tracks the module closely; near 0, it barely belongs.

The modules are tested two ways, as in `03_Set_Contrasts` for gene sets:

- Eigengene: the module as one variable. limma fits the eigengenes with the protein design (six
  cell means, blood slope, subject block) and asks whether the summary moved.
- fry: the module as a set of proteins. It asks whether the members moved together away from
  zero, with people as the unit of evidence, on limpa's complete matrix.

The two weight members differently. The eigengene weights each member by its loading on the
first component, so central proteins count most; fry weights every member equally. They can
therefore disagree on the same module.

Membership against significance asks whether a module's central proteins are the ones that moved:
within each module, the Spearman correlation between each member's kME and its protein-level
moderated t in the contrast. Positive means the most central proteins moved up most. Members of a
module are correlated by construction, so this correlation describes the pattern; its p is not a
test.

## What the code does

| Step | Function | Detail |
|---|---|---|
| Eigengene fit | `limma::duplicateCorrelation()`, `lmFit(block = subject)`, `contrasts.fit()`, `eBayes(robust = TRUE)` | Design and contrasts from `design.rds`; the within-subject correlation is re-estimated on the eigengenes. |
| Eigengene effects | `limma::topTable(confint = TRUE)` | logFC with its 95% interval, moderated t, p, BH FDR over modules within contrast, Storey's q (`storey_q()`, which returns BH when fewer than 100 modules are tested). |
| fry | `limma::fry(block = subject, correlation = ...)` | Each non-grey module as a protein set on limpa's complete matrix, correlation re-estimated on that matrix; FDR over modules within contrast. |
| Membership | `cor.test(method = "spearman", exact = FALSE)`, `DescTools::CorCI()` | kME against the protein moderated t from `fit.rds`, per module and contrast, with a Fisher-z 95% interval. |

## How to read the outputs

- With few modules, BH has little to adjust over, but each module is still one of several tests
  per contrast. A module at nominal p < 0.05 in one contrast is one result among many.
- An eigengene logFC with an interval crossing zero says the data allow no change; read the
  interval, not only p.
- Membership rho near zero with a significant eigengene means the module moved as a whole without
  its centre leading. A large rho in a module that did not move says little.
- kME was computed on subject-centred data in the trunk; the protein t comes from the full model.
  The two are on different footings by design.

## Files

| | |
|---|---|
| Reads | `01_Trunk/06_Module_Building/c_data/modules.rds` (eigengenes, membership, kME, sample metadata), `02_HRvLR_groups/01_Design_Contrasts/c_data/design.rds`, `02_HRvLR_groups/02_Protein_Contrasts/c_data/fit.rds` (protein moderated t), `01_Trunk/03_Quantification/c_data/DAList_limpa.rds` (complete matrix), `00_Input/phenotype.csv` (heatmap tracks) |
| Writes | `c_data/module_contrasts.xlsx` (read by `Figures/F2_Programmes`), `b_reports/module_contrasts_report.pdf` |
| Run | `Rscript 02_HRvLR_groups/04_Module_Contrasts/a_script/module_contrasts.R` |

```
04_Module_Contrasts/
├── a_script/module_contrasts.R
├── b_reports/module_contrasts_report.pdf
└── c_data/module_contrasts.xlsx
```

## Workbook: `module_contrasts.xlsx`

`read_me` lists the sheets; `columns` defines every column; `input_manifest` and
`package_versions` close the workbook.

| Sheet | One row is | Key columns |
|---|---|---|
| `eigengene_tests` | a module in a contrast | `logFC`, `logFC_lower`, `logFC_upper`, `t`, `p`, `fdr`, `qvalue` |
| `fry_tests` | a module in a contrast | `n_proteins`, `direction`, `p`, `fdr`, `qvalue` |
| `membership_contrast` | a module in a contrast | `n_members`, `rho`, `rho_lower`, `rho_upper`, `p` |
| `correlation` | a matrix | `matrix` (eigengenes or limpa), `correlation` |

## Report: `module_contrasts_report.pdf`

1. Every module's eigengene across the 45 biopsies (not centred), grouped by group and time point,
   with fCSA mixed, mCSA, leg press 1RM and comp_hypertrophy above (blank at T3). Rows follow
   module size.
2. Eigengene logFC in the nine contrasts beside a p/FDR status map.
3. Eigengene moderated t per module and contrast; label gives p where below 0.05, black border
   marks BH < 0.05.
4. fry signed -log10 p per module and contrast, same marks.
5. Membership rho per module and contrast.
6. kME against moderated t in `Training_Interaction`, one facet per module, with a
   least-squares line for reference.

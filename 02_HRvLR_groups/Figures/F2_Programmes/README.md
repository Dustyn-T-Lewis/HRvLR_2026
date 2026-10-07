# F2_Programmes

Parked draft: the script runs and its outputs are rebuilt, but the figure is not final.

Figure 2: which co-expression programmes move with training and the bout, and do the groups differ
in them?

## The idea

A programme here is a WGCNA module: proteins that rise and fall together across biopsies
(`01_Processing/06_Module_Building`). The figure says what each module is, how it behaves in the nine
contrasts, how the gene-set calls compare between fgsea and fry, and how one module looks person
by person.

| Panel | Draws | Rule that picks the rows |
|---|---|---|
| a | Each module's size as a bar in its WGCNA colour, named by its top over-represented gene set (or "no set at FDR < 0.05"), with its three top members (highest kME). | every module, largest first |
| b | Each module's eigengene in the nine contrasts, same rows as a. Tile fill is the eigengene logFC; a filled dot marks limma's eigengene FDR < 0.05, an open dot p < 0.05 only; a square marks fry p < 0.05 on the module's proteins. Columns sit under Training, Acute, Between groups and Interaction strips. | every module and contrast |
| c | Gene sets called per contrast, stacked up and down, in three columns: fgsea at FDR < 0.05, fry at FDR < 0.05, fry at p < 0.05. The dashed line in the last column is 5% of the sets tested, the count chance gives at p < 0.05. | every set and contrast |
| d | Each person's eigengene at T1, T2 and T3 for one module, one facet per group, with the group mean as a thick line. The subtitle gives that module's logFC and p in both interaction contrasts. | the module with the smallest eigengene p in either interaction contrast |

The eigengene is the first principal component of a module's members, one number per biopsy; its
logFC is a difference on the eigengene's own standardised scale. kME is a protein's correlation
with its module's eigengene. fgsea ranks a set against the other proteins and treats each protein
as an independent draw; fry asks whether the set's proteins moved away from zero across people,
under the design and subject block. Panel c puts the two side by side because they answer
different questions and can disagree.

## What the code does

| Step | Function | Detail |
|---|---|---|
| Panel a | `panel_module_identity()`, `module_identity()` in `R/figure_panels.R` | Labels, top set and top members from `module_characterisation.xlsx`, sheet `labels`. A top set is named only at FDR < 0.05. |
| Panel b | `panel_module_contrasts()` | `eigengene_tests` and `fry_tests` from `module_contrasts.xlsx`; marks from `evidence_level()` in `R/figure_theme.R`. |
| Panel c | `panel_set_calls()` | `set_tests` from `set_tests.rds`; chance line at `n_tested * 0.05`. |
| Panel d | `panel_module_people()` | Eigengenes and metadata from `modules.rds`. |
| Assembly | `patchwork`, `save_figure()` | 300 x 175 mm, `cairo_pdf`. |
| Checks | `stopifnot()` | Every module in the labels has eigengene tests; panel b has one tile per module and contrast; panel c's fgsea and fry FDR counts equal `set_summary`; panel d covers every person. |

## How to read the outputs

- Panel d's module is chosen for having the smallest interaction p among all modules, so its
  difference is inflated by the choice. It shows how that one result spreads across people, not
  that it is real.
- Panel a's subtitle states the module count as fixed text; it does not follow 01_Processing if the
  module count changes.
- A module's top set comes from over-representation among its members; it names the module, it
  does not test it.

## Files

| | |
|---|---|
| Reads | `02_HRvLR_groups/01_Design_Contrasts/c_data/design_contrasts.xlsx` (`contrasts`, for column order), `02_HRvLR_groups/04_Module_Contrasts/c_data/module_contrasts.xlsx` (`eigengene_tests`, `fry_tests`), `02_HRvLR_groups/03_Set_Contrasts/c_data/set_tests.rds`, `02_HRvLR_groups/03_Set_Contrasts/c_data/set_contrasts.xlsx` (`set_summary`, for the check), `01_Processing/06_Module_Building/c_data/modules.rds`, `01_Processing/07_Module_Characterisation/c_data/module_characterisation.xlsx` (`labels`) |
| Writes | `b_reports/f2_programmes.pdf`, `c_data/f2_programmes.xlsx` |
| Run | `Rscript 02_HRvLR_groups/Figures/F2_Programmes/a_script/f2_programmes.R` |

```
F2_Programmes/
├── a_script/f2_programmes.R
├── b_reports/f2_programmes.pdf
└── c_data/f2_programmes.xlsx
```

## Workbook: `f2_programmes.xlsx`

`read_me` lists the sheets; `columns` defines every column; `input_manifest` and
`package_versions` close the workbook.

| Sheet | One row is | Key columns |
|---|---|---|
| `panel_sources` | a panel | `panel`, `shows`, `file`, `sheet`, `rule` |
| `panel_a` | a module | `n_proteins`, `top_set`, `top_set_fdr`, `top_members`, `string_ratio` |
| `panel_b` | a module in a contrast | `logFC`, `p`, `fdr`, `p_fry`, `fdr_fry`, `family`, `column` |
| `panel_c` | a call and contrast | `call`, `family`, `row`, `contrast`, `n` |
| `panel_d` | a biopsy of the drawn module | `subject`, `arm`, `timepoint`, `value`, `module` |

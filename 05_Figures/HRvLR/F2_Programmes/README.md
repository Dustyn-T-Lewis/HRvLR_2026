# F2_Programmes

Parked draft: the script runs and its outputs are rebuilt, but the figure is not final.

Figure 2: which co-expression programmes (WGCNA modules) move with training and the bout, and do the groups differ in them?

| Panel | Draws | Rule that picks the rows |
|---|---|---|
| a | Each module's size as a bar in its WGCNA colour, named by its top over-represented gene set (or "no set at FDR < 0.05"), with its three top members (highest kME). | every module, largest first |
| b | Each module's eigengene in the nine contrasts, same rows as a. Tile fill is the eigengene logFC; a filled dot marks limma's eigengene FDR < 0.05, an open dot p < 0.05 only; a square marks fry p < 0.05 on the module's proteins. Columns sit under Training, Acute, Between groups and Interaction strips. | every module and contrast |
| c | Gene sets called per contrast, stacked up and down, in three columns: fgsea at FDR < 0.05, fry at FDR < 0.05, fry at p < 0.05. The dashed line in the last column is 5% of the sets tested, the count chance gives at p < 0.05. | every set and contrast |
| d | Each person's eigengene at T1, T2 and T3 for one module, one facet per group, with the group mean as a thick line. The subtitle gives that module's logFC and p in both interaction contrasts. | the module with the smallest eigengene p in either interaction contrast |

| | |
|---|---|
| Reads | `03_Group_comparisons/HRvLR/c_data/01_design_contrasts.xlsx` (`contrasts`, for column order), `03_Group_comparisons/HRvLR/c_data/04_module_contrasts.xlsx` (`eigengene_tests`, `fry_tests`), `03_Group_comparisons/HRvLR/c_data/03_set_tests.rds`, `03_Group_comparisons/HRvLR/c_data/03_set_contrasts.xlsx` (`set_summary`, for the check), `01_Processing/06_Module_Building/c_data/modules.rds`, `01_Processing/07_Module_Characterisation/c_data/module_characterisation.xlsx` (`labels`) |
| Writes | `b_reports/f2_programmes.pdf`, `c_data/f2_programmes.xlsx` |
| Run | `Rscript 05_Figures/HRvLR/F2_Programmes/a_script/f2_programmes.R` |

## Read with care

- Panel d's module has the smallest interaction p of all modules, so its difference is inflated by the choice; it shows how one result spreads across people.
- A module's top set comes from over-representation among its members; it names the module and tests nothing.
- Panel c sets fgsea (ranks a set against the other proteins) beside fry (do its proteins move across people); they answer different questions and can disagree.
- The script stops if a plotted count differs from `set_summary` or a module lacks eigengene tests.

## Sheets

| Sheet | One row is | Key columns |
|---|---|---|
| `panel_sources` | a panel | `panel`, `shows`, `file`, `sheet`, `rule` |
| `panel_a` | a module | `n_proteins`, `top_set`, `top_set_fdr`, `top_members`, `string_ratio` |
| `panel_b` | a module in a contrast | `logFC`, `p`, `fdr`, `p_fry`, `fdr_fry`, `family`, `column` |
| `panel_c` | a call and contrast | `call`, `family`, `row`, `contrast`, `n` |
| `panel_d` | a biopsy of the drawn module | `subject`, `arm`, `timepoint`, `value`, `module` |

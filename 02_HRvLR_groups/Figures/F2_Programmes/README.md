# F2_Programmes

Figure 2: which co-expression programmes move with training and the bout, and do HR and LR differ
in them?

| | |
|---|---|
| Reads | `02_HRvLR_groups/01_Design_Contrasts/c_data/design_contrasts.xlsx` (`contrasts`), `02_HRvLR_groups/04_Module_Contrasts/c_data/module_contrasts.xlsx` (`eigengene_tests`, `fry_tests`), `02_HRvLR_groups/03_Set_Contrasts/c_data/set_tests.rds`, `02_HRvLR_groups/03_Set_Contrasts/c_data/set_contrasts.xlsx` (`set_summary`), `01_Trunk/06_Module_Building/c_data/modules.rds`, `01_Trunk/07_Module_Characterisation/c_data/module_characterisation.xlsx` (`labels`) |
| Writes | `b_reports/f2_programmes.pdf` (300 x 175 mm), `c_data/f2_programmes.xlsx` |
| Run | `Rscript 02_HRvLR_groups/Figures/F2_Programmes/a_script/f2_programmes.R` |

The script reads saved outputs only. Before writing, it checks that every module and contrast is
drawn once and that the fgsea and fry FDR counts in panel c equal `set_summary`.

## Panels

| Panel | Shows | Rule |
|---|---|---|
| a | Each of the 15 trunk modules: size in its WGCNA colour, its top over-represented set (or "no set at FDR < 0.05"), and its three strongest hub proteins. | every module, largest first |
| b | Each module's eigengene on the nine contrasts, in the same rows as a. Fill is the eigengene logFC; a dot marks limma's eigengene p < 0.05 (filled at FDR < 0.05); a square marks fry p < 0.05 on the module's proteins. | every module and contrast |
| c | Gene sets called per contrast: fgsea at FDR < 0.05, fry at FDR < 0.05, and fry at p < 0.05, up and down. The dashed line is the fry count chance gives at p < 0.05 (5% of 1,378). | every set and contrast |
| d | Each person's eigengene at T1, T2 and T3 for the module with the smallest eigengene p in either interaction contrast, one facet per group, with group means. | smallest interaction p |

fgsea and fry ask different questions. fgsea ranks a set against the other proteins and counts
each protein as an independent draw; fry asks whether the set's proteins moved away from zero
across these people, under the same design and subject block. fgsea calls up to 299 sets in a
contrast and fry none at FDR < 0.05 in any, so fgsea counts are read as a ranking, not as findings.

No module eigengene reaches FDR < 0.05 on any contrast (smallest 0.092). The acute bout in HR
carries most of the nominal module signal; fry on the purple module's proteins in `Acute_HR` is
the one module test at FDR < 0.05 (0.043). The lead interaction module is salmon (cytoskeleton
organisation), acute interaction p = 0.037, FDR 0.48: a single nominal result, drawn so its
spread across people can be seen.

## Workbook

| Sheet | Holds |
|---|---|
| `panel_sources` | Per panel: what it shows, the file and sheet it reads, and the rule that picks its rows. |
| `panel_a` | Module size, top set and its FDR, three hub proteins, STRING enrichment ratio. |
| `panel_b` | Eigengene logFC, p and FDR, and fry p and FDR, per module and contrast. |
| `panel_c` | Sets called up and down per contrast by fgsea (FDR), fry (FDR) and fry (p). |
| `panel_d` | Each person's eigengene for the lead interaction module. |

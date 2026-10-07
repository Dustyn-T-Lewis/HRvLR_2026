# 01 · Trunk

Everything every arc reads, built without any group label.

| Substage | Does | Writes for the arcs |
|---|---|---|
| [`01_Protein_Filtering`](01_Protein_Filtering/README.md) | removes contaminants by identity, the three runs that fail on depth, and proteins seen too rarely | `DAList_filtered.rds` |
| [`02_Normalisation`](02_Normalisation/README.md) | cyclic loess | `DAList_normalized.rds` |
| [`03_Quantification`](03_Quantification/README.md) | limpa: a value and standard error for every protein in every sample; a blood-corrected copy | `proteins.rds`, `DAList_limpa.rds`, `DAList_corrected.rds` |
| [`04_Gene_Sets`](04_Gene_Sets/README.md) | freezes the MSigDB and GO-slim sets and maps proteins to genes | `gene_sets.rds` |
| [`05_Set_Scores`](05_Set_Scores/README.md) | one singscore per set and sample | `singscore.rds` |
| [`06_Module_Building`](06_Module_Building/README.md) | WGCNA co-expression modules and their eigengenes | `modules.rds` |
| [`07_Module_Characterisation`](07_Module_Characterisation/README.md) | names each module by its enriched sets, hubs and STRING edges | workbook only |
| [`Figures/S01_QC`](Figures/S01_QC/README.md) | quality-control supplement | figure |

Two choices keep the trunk free of labels. The blood slope is removed against time point only
(`removeBatchEffect`, design `~ 0 + timepoint`), and modules are defined on abundance centred within
each subject. Keeping the three failed runs was tried on 2026-10-06 and reversed: they lowered
limpa's detection-curve slope from 0.67 to 0.59, and limma's sample weights did not single them
out.

# F03 · Pathways

| | |
|---|---|
| Reads | `03_Pathway_Enrichment/01_run_fgsea_and_fry/c_data/01_run_fgsea_and_fry.xlsx`, `03_Pathway_Enrichment/03_enrich_scatter_fgsea/c_data/03_enrich_scatter_fgsea.xlsx` |
| Writes | `b_reports/F03.pdf`, `b_reports/F03.png`, one PDF and PNG per panel in `b_reports/panels/`, `c_data/F03_data.xlsx` |
| Run | `Rscript 05_Figures/F03/a_script/F03.R` (each panel also runs alone) |

- A: sets at FDR < 0.05 per contrast: fgsea, fgsea after collapse, fry.
- B: HR against LR training NES per set.
- C: ten strongest collapse survivors on Training_Interaction.

| File | Figure | Panel |
|---|---|---|
| `b_reports/F03.pdf`, `.png` | F03 | all |
| `b_reports/panels/A_set_counts.pdf`, `.png` | F03 | A |
| `b_reports/panels/B_nes_training.pdf`, `.png` | F03 | B |
| `b_reports/panels/C_top_sets.pdf`, `.png` | F03 | C |
| `c_data/F03_data.xlsx`, sheet `A_set_counts` | F03 | A |
| `c_data/F03_data.xlsx`, sheet `B_nes_training` | F03 | B |
| `c_data/F03_data.xlsx`, sheet `C_top_sets` | F03 | C |

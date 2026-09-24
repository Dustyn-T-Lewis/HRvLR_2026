# F04 · Networks

| | |
|---|---|
| Reads | `04_Network/01_build_modules/c_data/01_build_modules.xlsx`, `02_characterise_modules/c_data/02_characterise_modules.xlsx`, `03_preserve_modules/c_data/03_preserve_modules.xlsx`, `04_test_modules/c_data/04_test_modules.xlsx` |
| Writes | `b_reports/F04.pdf`, `b_reports/F04.png`, one PDF and PNG per panel in `b_reports/panels/`, `c_data/F04_data.xlsx` |
| Run | `Rscript 05_Figures/F04/a_script/F04.R` (each panel also runs alone) |

- A: proteins per module and eigengene subject ICC.
- B: STRING observed over expected edges, with each module's top set at FDR < 0.05.
- C: Zsummary of each arm's modules in the other arm.
- D: moderated t of each eigengene on the nine contrasts.

| File | Figure | Panel |
|---|---|---|
| `b_reports/F04.pdf`, `.png` | F04 | all |
| `b_reports/panels/A_modules.pdf`, `.png` | F04 | A |
| `b_reports/panels/B_string.pdf`, `.png` | F04 | B |
| `b_reports/panels/C_preservation.pdf`, `.png` | F04 | C |
| `b_reports/panels/D_eigengene_tests.pdf`, `.png` | F04 | D |
| `c_data/F04_data.xlsx`, sheet `A_modules` | F04 | A |
| `c_data/F04_data.xlsx`, sheet `B_string` | F04 | B |
| `c_data/F04_data.xlsx`, sheet `C_preservation` | F04 | C |
| `c_data/F04_data.xlsx`, sheet `D_eigengene_tests` | F04 | D |

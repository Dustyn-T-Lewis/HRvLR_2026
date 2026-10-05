# F03_Continuum

Figure 3: growth and strength as continuous traits, with no HR/LR labels.

| | |
|---|---|
| Reads | `02_Proteins/04_Associate/c_data/04_associate.xlsx` (`chance`, `interaction`, `change_score`), `03_Pathways/05_Associate/c_data/05_associate.xlsx` (`interaction`), `04_Modules/04_Associate/c_data/04_associate.xlsx` (`interaction`), `05_Integration/04_Continuous/c_data/04_continuous.xlsx` (`pathway_chance`, `composition`), `05_Integration/05_Clustering/c_data/05_clustering.xlsx` (`summary`), `04_Modules/01_Build/c_data/modules.rds`, `01_Preprocess/03_Quantify/c_data/DAList_corrected.rds`, `03_Pathways/01_Gene_Sets/c_data/gene_sets.rds`, `00_Input/phenotype.csv` |
| Writes | `b_reports/f03_continuum.pdf`, `b_reports/f03_continuum.png` (300 x 190 mm, 600 dpi), `c_data/f03_continuum.xlsx` |
| Run | `Rscript 06_Figures/F03_Continuum/a_script/f03_continuum.R` |

The script reads saved outputs only and checks its plotted numbers against them before writing.

## Panels

| Panel | Shows | Reading |
|---|---|---|
| a | Every continuous screen's count as a ratio to its shuffle 95th percentile (200 phenotype-row shuffles); proteins count p < 0.05, pathways count sets at FDR < 0.05 | 1 of 72 screens beats chance: T1 pathway levels against later leg-extension gain (p = 0.03) |
| b | The 20 proteins at BH whose training change depends on growth x leg-press gain, before and after adjusting for fibre composition | all sarcomere proteins (MYOT, MYH7, MYH1, MYBPC2, MYL3), all negative; 13 keep FDR < 0.05 after adjustment; TPM3, TTN, TNNI1, TNNT1, NEB, MYOM1 and PDLIM5 do not |
| c | Black module training change against leg-press gain, coloured by growth; the people whose removal weakens the protein results most are named, with how many proteins each | LR S30 for 14 proteins and HR S10 for 6; module FDR 0.0045 |
| d | The 23 pathways at BH in the same interaction | sarcomere organisation, myofibril assembly and striated muscle contraction lead; 22 of 23 negative |
| e | ENO2 and UQCR11 at T1 against leg-press 1RM at T1 | rho +0.90 and -0.90, FDR 0.026; the strongest hit beats the shuffles (familywise p = 0.015), but the screen's count does not |
| f | Clustering of people per input against a max-statistic null, and the best agreement of pathway-based groups with HR/LR | people cluster more tightly than chance for most inputs; no grouping matches the arms |

Every protein in b keeps p < 0.05 when any one person is left out (largest p 0.008). The result
still leans on two people, named in c, and on fibre composition for 7 of 20 proteins.

## Workbook

| Sheet | Holds |
|---|---|
| `screens` | Panel a: per screen, observed count, shuffle 95th percentile, ratio and permutation p. |
| `interaction_proteins` | Panel b: FDR before and after composition, leave-one-out maximum p and the person giving it. |
| `black_change` | Panel c: each person's module change, leg-press gain, growth and driver count. |
| `interaction_pathways` | Panel d: pathways at BH. |
| `baseline_strength` | Panel e: ENO2 and UQCR11 at T1 and leg-press 1RM at T1. |
| `clustering` | Panel f: clustering statistics and their nulls. |

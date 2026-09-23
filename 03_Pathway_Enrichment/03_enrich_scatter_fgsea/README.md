# 03_enrich_scatter_fgsea

Each set's HR NES against its LR NES, once for training and once for the acute bout.

| | |
|---|---|
| **Reads** | `set_tests.rds` |
| **Writes** | 4 composites, `nes_scatter.csv`, `03_enrich_scatter_fgsea.xlsx` |

| Pair | Population | Sets | rho | Significant | Discordant |
|---|---|---:|---:|---:|---:|
| training | all collections | 1,369 | 0.34 | 176 | 15 |
| training | collapse survivors | 38 | 0.66 | 38 | 6 |
| training | Hallmark and GO Slim | 90 | 0.28 | 11 | 3 |
| acute | all collections | 1,370 | 0.35 | 287 | 53 |
| acute | collapse survivors | 95 | 0.58 | 95 | 19 |
| acute | Hallmark and GO Slim | 90 | 0.49 | 43 | 10 |

Only sets scored in both contrasts are compared; a set near the size floor can drop out of one
ranking. A set is discordant when its two NES differ in sign.

Files are numbered by pair, `01_*` training and `02_*` acute. `*_all_*`: all sets, collapse
survivors, discordant sets. `*_curated_*`: Hallmark and GO Slim, then each concordant quadrant
rescaled so every significant set is named. These two collections are used because their sets do
not nest.

# 03_enrich_scatter_fgsea

Each set's HR NES against its LR NES, once for training and once for the acute bout.

| | |
|---|---|
| **Reads** | `set_tests.rds` |
| **Writes** | 4 composites, `nes_scatter.csv`, `03_enrich_scatter_fgsea.xlsx` |

| Pair | Population | Sets | rho | Significant | Discordant |
|---|---|---:|---:|---:|---:|
| training | all collections | 1,111 | 0.34 | 157 | 13 |
| training | collapse survivors | 32 | 0.63 | 32 | 6 |
| training | Hallmark and GO Slim | 58 | 0.14 | 7 | 3 |
| acute | all collections | 1,087 | 0.35 | 194 | 37 |
| acute | collapse survivors | 68 | 0.57 | 68 | 14 |
| acute | Hallmark and GO Slim | 55 | 0.46 | 24 | 6 |

Only sets scored in both contrasts are compared; a set near the size floor can drop out of one
ranking. A set is discordant when its two NES differ in sign.

`01_*` composites: all sets, collapse survivors, discordant sets. `02_*`: Hallmark and GO Slim,
then each concordant quadrant rescaled so every significant set is named. These two collections
are used because their sets do not nest.

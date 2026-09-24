# 03_preserve_modules

Whether HR and LR share the same co-expression structure.

| | |
|---|---|
| Reads | `01_Preprocess/03_Imputation/c_data/DAList_imputed.rds`, `01_build_modules/c_data/modules.rds` |
| Writes | `module_preservation.rds`, `03_preserve_modules.xlsx`, 1 figure |

Modules are built inside each arm with the settings of `01_build_modules`: subject-centred
abundance, signed network, bicor. HR gives 21 modules at power 16 (the WGCNA FAQ fallback for 20
samples, since no power cleared 0.85) and LR 17 at power 9. `WGCNA::modulePreservation()` then
tests each arm's modules in the other arm, 200 permutations.

Zsummary above 10 is strong preservation, 2 to 10 moderate, below 2 none (Langfelder et al.
2011). Median rank orders modules against each other and does not depend on size. Each arm module
is labelled with its best-overlapping full-cohort module.

| Direction | Strong | Moderate | None |
|---|---:|---:|---:|
| HR modules in LR | 3 | 17 | 1 |
| LR modules in HR | 4 | 11 | 2 |

Turquoise is strongly preserved both ways (Zsummary 28 and 31), as is the muscle contraction module
(HR purple 15.3 in LR; LR greenyellow, which maps to full-cohort purple, 16.2 in HR). Three modules of
39 to 54 proteins fall below 2: HR lightyellow, and LR lightcyan and tan.

Each arm has 20 to 24 samples, so these are small networks; read the result as a comparison of
structure, not as module discovery.

# 01_build_modules

WGCNA modules and their eigengenes.

| | |
|---|---|
| **Reads** | `01_Preprocess/03_Imputation/c_data/DAList_imputed.rds` |
| **Writes** | `modules.rds`, `01_build_modules.xlsx`, 3 figures |

Modules are **defined** on abundance centred within subject and **scored** on raw abundance. On
raw abundance subject identity drives the leading components, so modules built there would encode
who a biopsy came from. Centring removes that and leaves how proteins move together inside a
person; scoring on raw abundance keeps the between-arm differences testable.

Signed network, biweight midcorrelation (`maxPOutliers = 0.05`), `deepSplit = 2`,
`minModuleSize = 30`, `mergeCutHeight = 0.15`. The soft power is the lowest whose signed
scale-free fit clears 0.85: power 8, R2 0.874, mean connectivity 28. bicor warns about proteins
with zero MAD, which imputation can produce; `pearsonFallback = "individual"` handles them.

| Module | Proteins | Subject ICC |
|---|---:|---:|
| turquoise | 359 | 0.00 |
| blue | 174 | 0.00 |
| brown | 164 | 0.11 |
| yellow | 159 | 0.14 |
| green | 156 | 0.00 |
| red | 143 | 0.20 |
| black | 134 | 0.00 |
| pink | 81 | 0.37 |
| magenta | 76 | 0.12 |
| purple | 68 | 0.24 |
| greenyellow | 63 | 0.00 |
| tan | 58 | 0.37 |

Figures: soft-threshold curves, module sizes with subject ICC, eigengene trajectories by arm.

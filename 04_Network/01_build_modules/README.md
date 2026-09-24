# 01_build_modules

WGCNA modules and their eigengenes.

| | |
|---|---|
| Reads | `01_Preprocess/03_Imputation/c_data/DAList_imputed.rds` |
| Writes | `modules.rds`, `01_build_modules.xlsx`, 3 figures |

Modules are defined on abundance centred within subject and scored on raw abundance. On
raw abundance subject identity drives the leading components, so modules built there would encode
who a biopsy came from. Centring leaves how proteins move together inside a person; scoring on raw
abundance keeps between-arm differences testable. HR_S28 has one biopsy after outlier removal and
centres to zeros, so 44 samples define the modules and all 45 are scored.

Signed network, biweight midcorrelation (`maxPOutliers = 0.05`), `deepSplit = 2`,
`minModuleSize = 30`, `mergeCutHeight = 0.15`. The soft power is the lowest whose signed
scale-free fit clears 0.85: power 8, R2 0.866, mean connectivity 28. bicor warns about proteins
with zero MAD, which imputation can produce; `pearsonFallback = "individual"` handles them.

| Module | Proteins | Subject ICC |
|---|---:|---:|
| turquoise | 363 | 0.00 |
| blue | 182 | 0.17 |
| brown | 170 | 0.00 |
| yellow | 144 | 0.19 |
| green | 134 | 0.00 |
| red | 133 | 0.00 |
| black | 102 | 0.00 |
| pink | 94 | 0.00 |
| magenta | 93 | 0.06 |
| purple | 81 | 0.37 |
| greenyellow | 63 | 0.15 |
| tan | 56 | 0.35 |

Figures: soft-threshold curves, module sizes with subject ICC, eigengene trajectories by arm.

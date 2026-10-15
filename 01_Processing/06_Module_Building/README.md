# 06_Module_Building

Finds co-expression modules with WGCNA (signed, bicor) on person-centred abundance, then scores each
module's eigengene on the uncentred matrix so differences between people stay testable. No label or
phenotype enters construction.

## Inputs

| File | Holds |
|---|---|
| `01_Processing/03_Blood_Correction/c_data/proteins_corrected.rds` | limpa's complete matrix with blood slopes removed from observed values |

## Run

```sh
quarto render 01_Processing/06_Module_Building/a_script/module_building.qmd --output-dir ../b_reports
```

## Outputs

| File | Holds |
|---|---|
| `c_data/modules.rds` | `eigengenes` (module × biopsy), `membership` (module and kME per protein), `meta` (biopsy table), `module_summary`; read by every module-level step |
| `c_data/module_building.xlsx` | the same tables and the soft-power fit indices |
| `b_reports/module_building.html` | the rendered walkthrough |

## Read with care

- Soft power is fixed at 12, the WGCNA FAQ value for signed networks with more than 40 samples.
  The usual scale-free rule is met only at power 19, where the network is almost empty (mean
  connectivity 1.3, 73% of proteins in no module); a noise simulation of the same shape shows the
  fit rising with power on its own. The module QC compares powers 10, 14 and 19 with 12.
- Power 12 gives 13 modules of 63 to 389 proteins; 283 of 2,059 proteins fit no module. At powers
  10 and 14 the median module keeps 68% and 61% of its proteins in its best match; at 19, 30%.
- Module colours are labels only. They change if `01_Processing` is rebuilt.
- Eigengenes are scored on the uncentred matrix. Person ICCs run from 0 to 0.29; yellow's is
  exactly 0 (the mixed model's between-person variance is 0), so its variation is all within
  person.

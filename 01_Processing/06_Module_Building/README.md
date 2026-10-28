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
| `c_data/module_building.xlsx` | the same tables, the soft-power fit indices, the power QC, the settings screen and each module's leave-one-person-out stability |
| `b_reports/module_building.html` | the rendered walkthrough |

## Read with care

- The usual scale-free rule is met only at power 19, where the network is almost empty (mean
  connectivity 1.3, 73% of proteins in no module); a noise simulation of the same shape shows the
  fit rising with power on its own. Power 10 with modules merged at eigengene correlation 0.75
  (`mergeCutHeight = 0.25`) is used because, among 32 settings screened, it gives modules that
  best re-form when one person is left out (sheet `settings_screen`).
- No setting reaches the YvO rule (median Jaccard >= 0.5 with 75% of modules re-forming). Here the
  median is 0.52 and 7 of 10 modules re-form (`module_stability`); yellow and magenta do not.
  Finer splits (deepSplit 3, minimum size 20) re-form less often.
- Power 10 gives 10 modules of 70 to 340 proteins; 229 of 2,059 proteins fit no module. At powers
  8, 12 and 14 the median module keeps 77%, 70% and 70% of its proteins in its best match.
- Module colours are labels only. They change if `01_Processing` is rebuilt.
- Eigengenes are scored on the uncentred matrix. Person ICCs run from 0 to 0.27; green's is
  exactly 0 (the mixed model's between-person variance is 0), so its variation is all within
  person.

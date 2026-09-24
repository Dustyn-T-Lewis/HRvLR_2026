# 01_Design

Attaches the design formula and the nine contrasts to the DAList and records the within-subject
correlation.

| | |
|---|---|
| Reads | `01_Preprocess/02_Normalization/c_data/DAList_normalized.rds` |
| Writes | `c_data/design.rds`, `c_data/01_design.xlsx` |
| Run | `quarto render 02_Differential_Expression/01_Design/a_script/01_design.qmd --output-dir ../b_reports` |
| Cost | about 7 s |

`design.rds` holds the DAList with design and contrasts attached, the nine contrast strings, the
role table and the name of the floor contrast.

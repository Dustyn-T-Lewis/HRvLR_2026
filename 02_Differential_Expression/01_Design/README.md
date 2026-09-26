# 01_Design

Attaches the design formula and the nine contrasts to the DAList and records the within-subject
correlation.

| | |
|---|---|
| Reads | `01_Preprocess/02_Normalization/c_data/DAList_normalized.rds` |
| Writes | `c_data/design.rds`, `c_data/01_design.xlsx` |
| Run | `Rscript 02_Differential_Expression/01_Design/a_script/01_design.R` |
| Cost | about 7 s |

`design.rds` holds the DAList with design and contrasts attached, the nine contrast strings, the
role table and the name of the floor contrast. The script stops if the design is rank deficient or
a column name fails `make.names()`, since `add_contrasts()` parses the contrast strings as R code.

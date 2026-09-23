# 04 · Network

Finds proteins that move together in this data, whether or not a curated database groups them,
and asks the same questions of those groups that `02` and `03` asked of proteins and sets.

| Step | Runs | Writes |
|---|---|---|
| [`01_build_modules`](01_build_modules/README.md) | WGCNA on subject-centred abundance | `modules.rds`, 3 figures |
| [`02_characterise_modules`](02_characterise_modules/README.md) | ORA against the 03 sets, hubs, STRING | `module_characterisation.rds`, 2 figures |
| [`03_classify_and_associate_modules`](03_classify_and_associate_modules/README.md) | nine contrasts, eight tasks, ten phenotypes | `module_results.rds`, 5 figures |

```sh
for s in 01_build_modules 02_characterise_modules 03_classify_and_associate_modules; do
  Rscript 04_Network/$s/a_script/$s.R
done
```

`02_characterise_modules` needs the STRING files in `00_Input/downloads/`; `00_Input/README.md`
has the command.

## Results

Twelve modules at soft power 8; 265 of 1,900 proteins stay unassigned. Several modules name a
clear biology: pink is striated muscle contraction (hubs TTN, NEB, TPM1), tan and yellow are
respiratory chain and aerobic respiration, magenta and purple are translation, greenyellow is
extracellular matrix. All twelve share more STRING edges than their degrees predict (1.1 to 10.1
times; FDR < 0.05).

No module clears BH in any contrast, classification task or phenotype association. Nominal hits
sit at chance: 5 of 108 contrast tests against 5.4 expected, 6 of 96 task tests against 4.8, and
10 of 360 association tests against 18.

# 04 · Network

Finds proteins that move together in this data, whether or not a curated database groups them,
asks whether the two arms share that structure, and puts the groups through the tests `02` and
`03` applied to proteins and sets.

| Step | Runs | Writes |
|---|---|---|
| [`01_build_modules`](01_build_modules/README.md) | WGCNA on subject-centred abundance | `modules.rds`, 3 figures |
| [`02_characterise_modules`](02_characterise_modules/README.md) | ORA against the 03 sets, hubs, STRING enrichment, hub networks | `module_characterisation.rds`, 3 figures |
| [`03_preserve_modules`](03_preserve_modules/README.md) | modules built in each arm, preservation in the other | `module_preservation.rds`, 1 figure |
| [`04_test_modules`](04_test_modules/README.md) | nine contrasts on eigengenes and with fry, membership against significance | `module_tests.rds`, 5 figures |
| [`05_classify_and_associate_modules`](05_classify_and_associate_modules/README.md) | eight tasks, ten phenotypes | `module_results.rds`, 5 figures |

```sh
for s in 01_build_modules 02_characterise_modules 03_preserve_modules 04_test_modules \
         05_classify_and_associate_modules; do
  Rscript 04_Network/$s/a_script/$s.R
done
```

About five minutes, three of them in `03_preserve_modules`. `02_characterise_modules` needs the
STRING files in `00_Input/downloads/`; `00_Input/README.md` has the command.

## Results

Twelve modules at soft power 8; 285 of 1,900 proteins stay unassigned. Eight have a top set at
FDR < 0.05: purple is striated muscle contraction (hubs TTN, NEB, TPM1), tan respiratory electron
transport, blue oxidative phosphorylation, brown glycolysis, red mRNA metabolism, greenyellow rRNA
processing, magenta translation initiation, pink epithelial-mesenchymal transition. All twelve
share more STRING edges than their degrees predict (1.1 to 10.1 times; FDR < 0.05).

The two arms share their co-expression structure. Of 21 modules built in HR, 20 are at least
moderately preserved in LR (Zsummary above 2), and 15 of 17 LR modules in HR. The muscle
contraction module is strongly preserved in both directions (Zsummary 15.3 and 16.2).

No eigengene clears BH in any contrast, task or phenotype. fry, testing each module as a protein
set, finds one: red rises in Acute_LR (FDR 0.045 over the twelve modules). Nominal counts sit near
chance: 5 of 108 eigengene contrast tests against 5.4 expected, 7 of 108 fry tests, 6 of 96 task
tests against 4.8, and 12 of 360 association tests against 18.

Within modules, membership tracks protein-level change in blue and green: blue members with high kME
move most in the training interaction (rho 0.32), green members in the opposite direction on
the acute interaction (rho −0.51). Members share a module, so these correlations are descriptive,
not tests.

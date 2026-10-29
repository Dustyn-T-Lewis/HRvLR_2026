# Modules: supplementary module pages

Per-module and cross-module pages from `02_HRvLR/08_Module_Detail`, drawn by
`a_script/modules_pages.R`. Modules are named as in F5 (`style/module_labels.R`).

| Page | Shows |
|---|---|
| `b_reports/<module>/hubs.pdf` | the module's 15 highest-kME proteins (most central at the top) in the training and acute contrasts, HR, LR and HR − LR: fill logFC (capped at ±1), size −log10 p, filled at p < 0.05, black ring at FDR < 0.05 |
| `b_reports/<module>/rewiring.pdf` | each member's kME within LR's 24 biopsies against within HR's; off the diagonal, a protein moves with the module in one group more than the other; the 8 largest shifts named |
| `b_reports/overview/coupling.pdf` | Spearman rho between every pair of modules' eigengene changes (training, acute) within HR and within LR, and HR − LR boxed at Fisher p < 0.05 |
| `b_reports/overview/per_person.pdf` | each person's eigengene change per module and window, people ordered by the growth score |
| `b_reports/overview/pls.pdf` | sPLS-DA (HR against LR) and sPLS-R (trait gains) on the eigengenes: held-out score (balanced, repeated 4-fold CV) against its permutation null, and the component-1 loadings |

## Run

```sh
Rscript run_all.R 04_Figures/Modules
```

## Outputs

| File | Holds |
|---|---|
| `b_reports/<module>/`, `b_reports/overview/` | the pages above |
| `c_data/Modules.xlsx` | the tables drawn, column dictionary, input md5 |

## Read with care

- Coupling: after the bout, HR's modules move as one block and LR's do not. Respiration (red)
  and cristae (magenta) change together in HR (rho 0.90) but not in LR (0.17); mRNA processing
  (blue), brown and cytoskeleton (green) change together in HR (0.88 to 0.93) and opposite to
  glycolysis (pink, -0.93 to -0.95), with no such pattern in LR. 7 of 90 pair-windows differ at
  Fisher p < 0.05, 6 of them after the bout. Each rho rests on 8 people.
- Rewiring: in the respiration module, NDUFA11 and the mitochondrially encoded MT-ND2 move with
  the module in HR (kME 0.73 and 0.82) and not in LR (-0.16 and 0.01); talin-2 (TLN2) does the same in
  the cristae module (0.78 against -0.21).
- PLS: no model beats its permutation null at p < 0.05. The closest is the training changes
  predicting the muscle-area gain (Q² 0.06, p 0.10); HR against LR is at chance from every window.
- No hub protein passes FDR < 0.05 in any contrast shown.

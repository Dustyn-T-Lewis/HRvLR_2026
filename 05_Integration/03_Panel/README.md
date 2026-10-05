# 03_Panel

Asks whether a multi-feature model can tell HR from LR in subjects it did not train on: an elastic
net and DIABLO. Exploratory: nothing here is validated outside this cohort.

| | |
|---|---|
| Reads | `01_Preprocess/03_Quantify/c_data/DAList_limpa.rds`, `04_Modules/01_Build/c_data/modules.rds`, `00_Input/phenotype.csv` |
| Writes | `c_data/03_panel.xlsx`, `b_reports/03_panel_figures.pdf` |
| Run | `Rscript 05_Integration/03_Panel/a_script/03_panel.R` (about an hour on eight cores) |

## Method

Ten panels: five inputs (T1, T2 and T3 levels, training change, acute change) by two feature sets
(the module eigengenes, and every protein of limpa's complete matrix, so nothing is filled in
inside a fold). The input is limpa's matrix with the blood slope removed against time
point only, never arm, so the labels cannot inform the correction; the eigengenes are rescored on
that matrix with the modules from `01_Build`. `panel_cv()` wraps `nestedcv::nestcv.glmnet()`:

| Setting | Value | Source |
|---|---|---|
| Model | penalised logistic regression, alpha 0.5 fixed | Zou & Hastie 2005 |
| Lambda | `lambda.1se` from an inner leave-one-out in each training fold | Varma & Simon 2006; Cawley & Talbot 2010 |
| Outer loop | stratified 5-fold, 10 repeats; score = mean out-of-fold AUC | Smith et al. 2014 |
| Scaling | glmnet standardises within each training fold | |
| Stability | share of the 50 repeat-by-fold models that kept each feature | Meinshausen & Bühlmann 2010 |
| Significance | 1,000 label permutations of the whole pipeline, p = (b + 1)/(m + 1) | Phipson & Smyth 2010 |

The outer loop is stratified k-fold, not leave-one-out. On pooled leave-one-out an empty model
scores AUC 0: each held-out HR subject leaves a training set tipped towards LR (Airola et al. 2011).
Stratified folds shrink that bias without removing it, so the permutation null centres near 0.39,
not 0.5. Judge each panel against its null, never against 0.5.

Any T3 panel inherits the blood confound. LR's T3 biopsies carry more blood, and the blood
adjustment fits its slope within group cells, so the between-arm blood difference at T3 stays in
the data; the T3 module panel's two modules track the blood index (salmon rho = -0.83).

The module panels carry a leak: WGCNA built the modules on all 16 subjects, including each fold's
held-out ones. The protein panels do not share it. Report against TRIPOD+AI (Collins et al. 2024):
internal validation only, no external cohort.

DIABLO (`mixOmics::block.splsda`) asks the same question with blocks that must agree: modules and
proteins, alone and with a phenotype block where the input has one (T1: baseline trait levels;
training: % trait change, composite and volume; T2: post-training levels). keepX is fixed at 5 per
block, one component, so nothing is tuned on the answer; 5 repeats of stratified 5-fold CV and 200
label permutations of the whole procedure. The training phenotype block holds fibre-area change,
the label's own ingredient, so that model's AUC sits near 1 by construction and is marked not
interpretable; what it shows is which modules and proteins the model keeps beside fibre change.

## Figures

Page one: each elastic-net panel's permutation null with the observed mean AUC and its p. Page two:
up to 15 features per panel by how often the model kept them. Page three: each DIABLO model's AUC
against its null.

## Workbook

| Sheet | Holds |
|---|---|
| `panels` | Per panel: subjects, features, mean and range of the 10 repeat AUCs, share of models that kept no feature, permutation null median and 95th percentile, permutation p. |
| `selection` | Per panel and feature: share of the 50 repeat-by-fold models that kept it. |
| `diablo` | Per DIABLO model: input, phenotype block, arm sizes, mean CV AUC, permutation null and p, whether the AUC is interpretable. |
| `diablo_selection` | Per DIABLO model and feature: share of the 25 CV fits that kept it. |

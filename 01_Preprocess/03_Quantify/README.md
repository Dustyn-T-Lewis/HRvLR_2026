# 03_Quantify

Models the missing values with limpa, so every protein has a value in every sample with a standard
error, and writes the matrices the later steps read.

| | |
|---|---|
| Reads | `02_Normalize/c_data/DAList_normalized.rds` |
| Writes | `c_data/proteins.rds`, `c_data/DAList_limpa.rds`, `c_data/DAList_corrected.rds`, `c_data/03_quantify.xlsx`, `b_reports/03_quantify_figures.pdf` |
| Run | `Rscript 01_Preprocess/03_Quantify/a_script/03_quantify.R` |

## Method

A protein goes undetected more often when it is faint: here 31% of values are missing in the
faintest quarter of proteins and 1% in the brightest. A missing value therefore carries
information, and both ignoring it and imputing it as if measured bias the result. limpa (Li,
Cobbold & Smyth 2025) fits one detection probability curve for the run with `dpc()`: the
probability that a value is detected, as a logistic function of its true intensity. Here the
slope is 0.67, a little below the 0.7 to 1 the authors find typical; they note it runs low when
values are noisy.

`dpcQuantByRow()` then gives every protein a value in every sample, with a posterior standard error
of about 0.01 log2 where observed and about 0.8 where missing. Our input is already summarised to
proteins, so each row is its own feature: the route the limpa vignette gives for protein-level
input ("Working with data already summarized at the protein-level"). The authors recommend
precursor-level input when it exists; without it there are no peptides to pool, so the gain over
other methods is smaller than in a precursor pipeline such as BFR_limpa.

Three objects come out:

| File | Holds | Read by |
|---|---|---|
| `proteins.rds` | limpa's EList (values, standard errors, observation counts) with the sample sheet | the protein fit (`dpcDE`), the leave-one-out reruns, the permutation |
| `DAList_limpa.rds` | the complete matrix | fry, the overview row order, the panel |
| `DAList_corrected.rds` | the complete matrix with each protein's blood slope removed (`limma::removeBatchEffect`, design = the six group cells) and the observed mask | singscore, WGCNA, the rank tests, correlations and ANCOVA (all on the complete matrix), and `observed_share` |

Blood is centred here: `blood` is the 01_Filter index (mean log2 HBB and HBA1 against the
sample median) minus its mean.

## Figures

One page: the detection probability curve.

## Workbook

| Sheet | Holds |
|---|---|
| `curve` | DPC intercept and slope on the logit scale, and the share of values missing. |
| `proteins` | Per protein: values observed, mean log2 abundance where observed and over the complete matrix, and the median standard error where observed and where missing. |
| `blood` | Per sample: blood index and its centred value. |
| `sample_correlation` | Pearson correlation between samples on the blood-corrected matrix, each protein centred on its mean so the correlation reflects how samples differ rather than the abundance order they share. |

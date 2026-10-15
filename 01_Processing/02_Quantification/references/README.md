# Contaminant list

`contaminants.csv` lists the 222 proteins removed as contaminants after quantification: 120
plasma, 41 immunoglobulin, 19 erythrocyte, 18 keratin, 12 complement, 7 globin and 5 leukocyte,
each with its UniProt accession, gene, class and reason.

The list is frozen from the earlier protein-level pipeline (`main`,
`01_Processing/01_Protein_Filtering`, sheet `contaminants_removed`). That step built it from a
curated blood list, the Human Protein Atlas (immunoglobulin class, "Secreted to blood", plasma
concentration) and a myonuclear-expression rescue that kept muscle proteins. Those rules are not
rerun here.

Proteins found only in the precursor data cannot be on the list. `quantification.qmd` prints, and
does not remove, any remaining keratin, immunoglobulin or haemoglobin gene for review.

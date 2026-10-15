# Module names for F5 and S5. WGCNA colours are labels only and change when the network is rebuilt,
# so a module is named by its leading set (its top over-represented set in
# 07_Module_Characterisation, when that set's FDR is below alpha_fdr) with the colour after it.
# Expects figure_style.R sourced.

module_identity <- function(labels, width = 34) {
  identity <- dplyr::if_else(
    labels$top_set_fdr < alpha_fdr, pathway_label(labels$top_set, width),
    sprintf("no set at FDR < %s", alpha_fdr)
  )
  rlang::set_names(paste0(identity, " (", labels$module, ")"), labels$module)
}

characterisation_input <-
  "01_Processing/07_Module_Characterisation/c_data/module_characterisation.xlsx"
characterisation_file <- here::here(characterisation_input)
# One row per module, largest first: the row order of every module panel.
read_module_labels <- function() {
  readxl::read_excel(characterisation_file, "labels") |>
    dplyr::arrange(dplyr::desc(n_proteins), module)
}

# Shared by the S6 panels: the elastic-net results of 03_Continuous/04_Elastic_Net and the lead
# they draw, picked by rule: the fit with the smallest permutation p (ties: the higher Q2); its
# pairing and trait are then shown for every feature set. Expects dplyr attached.
elastic_net_input <- "03_Continuous/04_Elastic_Net/c_data/elastic_net.rds"
elastic_net <- readRDS(here::here(elastic_net_input))
lead <- elastic_net$summary |>
  arrange(p_perm, desc(q2)) |>
  slice_head(n = 1) |>
  select(pairing, trait, window)
in_lead <- \(d) semi_join(d, lead, by = c("pairing", "trait"))
feature_set_names <- c(proteins = "Proteins", modules = "Modules", hallmark = "Hallmark sets")

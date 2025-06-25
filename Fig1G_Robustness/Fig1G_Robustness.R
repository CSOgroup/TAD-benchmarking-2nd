source("Fig1G_Robustness/utils_Robustness.R")

tads_input <- readr::read_tsv("Data/TADs.tsv") %>%
  dplyr::filter(meta.sub != ".")

res_robustness <- moc_multitool(tads_input)

res_robustness

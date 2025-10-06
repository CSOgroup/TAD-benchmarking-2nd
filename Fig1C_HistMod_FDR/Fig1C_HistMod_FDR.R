source("Fig1C_HistMod_FDR/utils_HistMod.R")

tads_input <- readr::read_tsv("Data/TADs.tsv") %>%
  dplyr::filter(meta.sub == ".")

res_histone <- analyze_histone_modification_all_tools(tads_input)

res_histone

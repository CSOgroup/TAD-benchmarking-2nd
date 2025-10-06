source("Fig1E_CTCF_orientation/utils_CTCF.R")

tads_input <- readr::read_tsv("Data/TADs.tsv") %>%
  dplyr::filter(meta.sub == ".")

res_ctcf <- analyze_ctcf_orientation_all_tools(tads_input)

res_ctcf

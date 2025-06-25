source("Fig1C_HistMod_FDR/utils_HistMod.R")

tads_input <- readr::read_tsv("Data/TADs.tsv") %>%
  dplyr::filter(meta.sub == ".")

gr_h27 <- gr_protein_gain("H3K27me3")
gr_h36 <- gr_protein_gain("H3K36me3")

res_histone <- analyze_histone_modification_all_tools(tads_input, gr_h27 = gr_h27, gr_h36 = gr_h36)

res_histone

source("Fig1D_StructProt_EnrichBoundaries/utils_StructProt.R")

tads_input <- readr::read_tsv("Data/TADs.tsv") %>%
  dplyr::filter(meta.sub == ".")

res_protein <- analyze_proteins_enrichment_all_tools(tads_input)

res_protein

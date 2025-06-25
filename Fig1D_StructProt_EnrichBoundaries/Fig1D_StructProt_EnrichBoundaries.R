source("Fig1D_StructProt_EnrichBoundaries/utils_StructProt.R")

tads_input <- readr::read_tsv("Data/TADs.tsv") %>%
  dplyr::filter(meta.sub == ".")

gr_ctcf <- gr_protein_gain("CTCF")
gr_rad21 <- gr_protein_gain("RAD21")
gr_smc3 <- gr_protein_gain("SMC3")

res_protein <- analyze_multiple_proteins_boundary_enrichment(tads_input, gr_ctcf = gr_ctcf,
                                                             gr_rad21 = gr_rad21, gr_smc3 = gr_smc3)

res_protein

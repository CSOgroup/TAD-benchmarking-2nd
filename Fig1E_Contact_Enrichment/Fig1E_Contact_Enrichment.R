source("Fig1E_Contact_Enrichment/utils_CE.R")

tads_input <- readr::read_tsv("Data/TADs.tsv") %>%
  dplyr::filter(meta.sub == ".")


# Attention:
# The HiC_matrix file is too large to include here. Please download the data directly from our website,
# TADShop: [https://tadshop.unil.ch/](https://tadshop.unil.ch/)

row_mx <- readr::read_tsv("HiC_matrix.txt",
                          col_names = F, col_types = list(readr::col_double())) %>%
  as.matrix()

enrichment_matrix <- calculate_contact_enrichment_matrix(row_mx)

res_ce <- calculate_contact_enrichment_for_multiple_tools(
  tads_input, enrichment_matrix = enrichment_matrix
)

res_ce

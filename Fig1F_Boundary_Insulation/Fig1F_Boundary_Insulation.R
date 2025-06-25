source("Fig1F_Boundary_Insulation/utils_Insulation.R")

tads_input <- readr::read_tsv("Data/TADs.tsv") %>%
  dplyr::filter(meta.sub == ".")


# Attention:
# The HiC_matrix file is too large to include here. Please download the data directly from our website,
# TADShop: [https://tadshop.unil.ch/](https://tadshop.unil.ch/)

row_mx <- readr::read_tsv("HiC_matrix.txt",
                          col_names = F, col_types = list(readr::col_double())) %>%
  as.matrix()

df_ISindex <- calculate_insulation_scores(row_mx)

res_is <- calculate_insulation_scores_for_multiple_tools(
  tads_input, gap_len = 3, bg_len = 30, df_ISindex = df_ISindex
)

res_is

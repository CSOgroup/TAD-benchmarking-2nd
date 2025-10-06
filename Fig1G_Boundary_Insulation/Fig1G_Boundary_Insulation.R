source("Fig1F_Boundary_Insulation/utils_Insulation.R")

tads_input <- readr::read_tsv("Data/TADs.tsv") %>%
  dplyr::filter(meta.sub == ".")


# Attention:
# The HiC_matrix file is too large to include here. Please download the data directly from our website,
# TADShop: [https://tadshop.unil.ch/](https://tadshop.unil.ch/)

chr_sels <- unique(tads_input$chr)
res_list <- vector("list", length = length(chr_sels))
names(res_list) <- chr_sels
for (chr_sel in chr_sels) {
  row_mx <- readr::read_tsv(sprintf("HiC_matrix_%s.txt", chr_sel),
                            col_names = F, col_types = list(readr::col_double())) %>%
    as.matrix()
  df_ISindex <- calculate_insulation_scores(row_mx)
  res_list[[chr_sel]] <- tads_input %>%
    dplyr::filter(chr %in% chr_sel) %>%
    calculate_insulation_scores_for_multiple_tools_genome(gap_len = 3, bg_len = 30, df_ISindex = df_ISindex)
}

res_is <- dplyr::bind_rows(res_list) %>%
  dplyr::group_by(meta.tool, meta.resol, meta.sub) %>%
  dplyr::summarise(is_sum = sum(is_sum), is_len = sum(is_len), .groups = "drop") %>%
  dplyr::mutate(is = is_sum / is_len) %>%
  dplyr::select(-is_sum, -is_len)

res_is

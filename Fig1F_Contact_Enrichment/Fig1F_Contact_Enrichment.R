source("Fig1E_Contact_Enrichment/utils_CE.R")

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
  enrichment_matrix <- calculate_contact_enrichment_matrix(row_mx)
  res_list[[chr_sel]] <- tads_input %>%
    dplyr::filter(chr %in% chr_sel) %>%
    calculate_contact_enrichment_for_multiple_tools_genome(enrichment_matrix = enrichment_matrix)
}

res_ce <- dplyr::bind_rows(res_list) %>%
  dplyr::group_by(meta.tool, meta.resol, meta.sub) %>%
  dplyr::summarise(contact_enrichment_sum = sum(contact_enrichment_sum), n = sum(n), .groups = "drop") %>%
  dplyr::mutate(contact_enrichment = contact_enrichment_sum / n) %>%
  dplyr::select(-contact_enrichment_sum, -n)

res_ce

source("utils/utils.R")

calculate_protein_enrichment_single_tool <- function(df_diffresol_one, gr_protein, chr_sel = "chr6", gr_flank = 10000) {
  chrSize <- chrsize_hg19()[[chr_sel]]
  smallbin <- 5000
  peak_bg_gap <- 10000

  tad_edges <- unique(c((df_diffresol_one$start - 1), df_diffresol_one$end))
  tad_edges_new <- tad_edges[
    abs(tad_edges - chrSize) > gr_flank & abs(tad_edges - 1) > gr_flank
  ]

  if (length(tad_edges_new) == 0) {
    return(
      df_diffresol_one %>%
        dplyr::select(meta.tool, meta.resol, meta.sub, chr) %>%
        dplyr::distinct() %>%
        dplyr::mutate(fc_over_bg = 0)
    )
  }
  gr_chr <- tibble::tibble(chr = chr_sel, start = 1, end = chrSize) %>%
    GenomicRanges::GRanges()
  gr_tads <- tibble::tibble(
    chr = chr_sel, start = tad_edges_new, end = tad_edges_new
  ) %>%
    GenomicRanges::GRanges() %>%
    magrittr::add(gr_flank)
  gr_bg <- GenomicRanges::setdiff(gr_chr, magrittr::add(gr_tads, peak_bg_gap))
  peak_val <- gr_tads %>%
    bin_trim(smallbin) %>%
    unlist() %>%
    GenomicRanges::countOverlaps(gr_protein)
  bg <- gr_bg %>%
    bin_trim(smallbin) %>%
    unlist() %>%
    GenomicRanges::countOverlaps(gr_protein)
  boundary_value <- mean(peak_val)
  background_value <- mean(bg)
  df_diffresol_one %>%
    dplyr::select(meta.tool, meta.resol, meta.sub, chr) %>%
    dplyr::distinct() %>%
    dplyr::mutate(
      boundary_value = boundary_value,
      background_value = background_value
    )
}


analyze_multiple_proteins_boundary_enrichment <- function(df_diffresol, gr_ctcf = NULL, gr_rad21 = NULL, gr_smc3 = NULL, chr_sel = "chr6", gr_flank = 10000) {
  df_diffresol_list <- df_diffresol %>%
    dplyr::mutate(grp = paste(meta.tool, meta.resol, meta.sub, chr, sep = "||")) %>%
    split(.$grp)

  proteins <- list(CTCF = gr_ctcf, RAD21 = gr_rad21, SMC3 = gr_smc3)

  res_all <- lapply(
    df_diffresol_list,
    function(df_i) {
      purrr::imap_dfr(
        proteins,
        ~calculate_protein_enrichment_single_tool(
          df_i, gr_protein = .x, chr_sel = chr_sel, gr_flank = gr_flank
        ) %>%
          dplyr::mutate(protein = .y)
      )
    }
  ) %>% dplyr::bind_rows()

  res_all
}

analyze_proteins_enrichment_all_tools <- function(tb_tads, chip_folder = "Data") {
  # chr_sels <- paste0("chr", c(1:22, "X"))
  chr_sels <- unique(tb_tads$chr)
  res_list <- vector("list", length = length(chr_sels))
  names(res_list) <- chr_sels
  for (chr_sel in chr_sels) {
    gr_ctcf <- gr_protein_gain("CTCF", chr_sel = chr_sel, chip_folder = chip_folder)
    gr_rad21 <- gr_protein_gain("RAD21", chr_sel = chr_sel, chip_folder = chip_folder)
    gr_smc3 <- gr_protein_gain("SMC3", chr_sel = chr_sel, chip_folder = chip_folder)

    res_list[[chr_sel]] <- tb_tads %>%
      dplyr::filter(chr %in% chr_sel) %>%
      analyze_multiple_proteins_boundary_enrichment(gr_ctcf = gr_ctcf, gr_rad21 = gr_rad21,
                                                    gr_smc3 = gr_smc3, chr_sel = chr_sel)
  }

  dplyr::bind_rows(res_list) %>%
    dplyr::group_by(meta.tool, meta.resol, meta.sub, protein) %>%
    dplyr::summarise(boundary_value = sum(boundary_value), background_value = sum(background_value), .groups = "drop") %>%
    dplyr::mutate(fc_over_bg = boundary_value / background_value) %>%
    dplyr::select(-boundary_value, -background_value)
}

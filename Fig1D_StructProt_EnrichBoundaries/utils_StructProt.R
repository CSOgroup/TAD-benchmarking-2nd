source("utils/utils.R")

calculate_protein_enrichment_at_boundary <- function(df_diffresol_one, gr_protein, chr_sel = "chr6") {
  chrSize <- chrsize_hg19()[[chr_sel]]
  gr_flank <- 10000
  smallbin <- 5000
  peak_bg_gap <- 10000

  meta.tool <- unique(df_diffresol_one$meta.tool)
  meta.resol <- unique(df_diffresol_one$meta.resol)

  tad_edges <- unique(c(df_diffresol_one[[2]] - 1, (df_diffresol_one[[3]])))
  tad_edges_new <- tad_edges[
    abs(tad_edges - chrSize) > gr_flank & abs(tad_edges - 1) > gr_flank
  ]

  if (length(tad_edges_new) == 0) {
    return(
      tibble::tibble(
        meta.tool = meta.tool,
        meta.resol = meta.resol,
        fc_over_bg = 0
      )
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
  fc_over_bg <- mean(peak_val) / mean(bg)
  tibble::tibble(
    meta.tool = meta.tool, meta.resol = meta.resol,
    fc_over_bg = fc_over_bg
  )
}


analyze_multiple_proteins_boundary_enrichment <- function(df_diffresol, gr_ctcf = NULL, gr_rad21 = NULL, gr_smc3 = NULL, chr_sel = "chr6") {
  df_diffresol_list <- df_diffresol %>%
    split(.$meta.tool)

  res_ctcf <- df_diffresol_list %>%
    lapply(calculate_protein_enrichment_at_boundary, gr_protein = gr_ctcf, chr_sel = chr_sel) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(protein = "CTCF")
  res_rad21 <- df_diffresol_list %>%
    lapply(calculate_protein_enrichment_at_boundary, gr_protein = gr_rad21, chr_sel = chr_sel) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(protein = "RAD21")
  res_smc3 <- df_diffresol_list %>%
    lapply(calculate_protein_enrichment_at_boundary, gr_protein = gr_smc3, chr_sel = chr_sel) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(protein = "SMC3")

  list(res_ctcf, res_rad21, res_smc3) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(resolution_kb = meta.resol / 1000) %>%
    dplyr::select(meta.tool, resolution_kb, protein, fc_over_bg)
}

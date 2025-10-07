source("utils/utils.R")

restrict_tads_to_ctcf <- function(gr_tads, gr_ctcf) {
  hits <- GenomicRanges::findOverlaps(gr_tads, gr_ctcf) %>%
    tibble::as_tibble()
  gr_intersect <- GenomicRanges::pintersect(gr_tads[hits$queryHits], gr_ctcf[hits$subjectHits])
  return(gr_intersect)
}


calculate_CTCF_orientation_at_boundary <- function(df_diffresol_one, gr_ctcf, chr_sel = "chr6", pwm, gr_flank = 10000) {
  chrSize <- chrsize_hg19()[[chr_sel]]
  # gr_flank <- 10000

  # tad_edges <- unique(c((df_diffresol_one$start - 1), df_diffresol_one$end))
  tad_edges <- unique(df_diffresol_one$tad_boundary)
  tad_edges_new <- tad_edges[
    abs(tad_edges - chrSize) > gr_flank & abs(tad_edges - 1) > gr_flank
  ]

  if (length(tad_edges_new) == 0) {
    return(
      tibble::tibble(
        chr = character(),
        tad_boundary = double(),
        ctcf_binding = double()
      )
    )
  }

  gr_tads <- tibble::tibble(
    chr = chr_sel, start = tad_edges_new, end = tad_edges_new, id = tad_edges_new
  ) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
    magrittr::add(gr_flank)
  gr_tads_ctcf <- restrict_tads_to_ctcf(gr_tads, gr_ctcf)

  if (length(gr_tads_ctcf) == 0) {
    res <- df_diffresol_one %>%
      dplyr::mutate(ctcf_binding = FALSE)
    return(res)
  }

  seqs <- Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19, gr_tads_ctcf)
  names(seqs) <- paste(GenomicRanges::start(gr_tads_ctcf), GenomicRanges::end(gr_tads_ctcf), sep = "_")

  hits_list <- lapply(seq_along(seqs), function(i) {
    TFBSTools::searchSeq(pwm, seqs[[i]], seqname = names(seqs)[i]) %>%
      as("GRanges") %>%
      tibble::as_tibble()
  }) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(s_e = as.character(seqnames)) %>%
    tidyr::separate(s_e, c("tad_ctcf_region_start", "tad_ctcf_region_end")) %>%
    dplyr::mutate(tad_ctcf_region_start = as.double(tad_ctcf_region_start), tad_ctcf_region_end = as.double(tad_ctcf_region_end)) %>%
    dplyr::mutate(ctct_binding_start = tad_ctcf_region_start + start - 1, ctct_binding_end = tad_ctcf_region_start + end - 1) %>%
    dplyr::select(tad_ctcf_region_start, tad_ctcf_region_end, ctct_binding_start, ctct_binding_end, strand, absScore, relScore, ID, siteSeqs)

  res <- gr_tads_ctcf %>%
    tibble::as_tibble() %>%
    dplyr::mutate(chr = as.character(seqnames)) %>%
    dplyr::select(chr, tad_ctcf_region_start = start, tad_ctcf_region_end = end, id) %>%
    dplyr::distinct() %>%
    dplyr::group_by(tad_ctcf_region_start, tad_ctcf_region_end, chr) %>%
    dplyr::summarise(id = paste(unique(id), collapse = ","), .groups = "drop") %>%
    dplyr::full_join(hits_list, by = c("tad_ctcf_region_start", "tad_ctcf_region_end")) %>%
    tidyr::separate_rows(id, sep = ",") %>%
    dplyr::select(chr, tad_boundary = id, strand, tad_ctcf_region_start, tad_ctcf_region_end, ctct_binding_start, ctct_binding_end, absScore, relScore, ID, siteSeqs) %>%
    dplyr::mutate(ctcf_binding = TRUE, tad_boundary = as.double(tad_boundary))

  tad_boundary_left <- setdiff(tad_edges, res$tad_boundary)
  if (length(tad_boundary_left) > 0) {
    res_left <- tibble::tibble(chr = chr_sel, tad_boundary = tad_boundary_left, ctcf_binding = FALSE)
    res <- dplyr::bind_rows(res, res_left)
  }

  res
}

analyze_ctcf_orientation_all_tools <- function(tb_all, chip_folder = "Data") {
  tb_tads <- tb_all %>%
    dplyr::select(chr, start, end) %>%
    dplyr::mutate(start = start - 1) %>%
    tidyr::pivot_longer(cols = start:end, names_to = "type", values_to = "tad_boundary") %>%
    dplyr::select(chr, tad_boundary) %>%
    dplyr::distinct()
  pfm <- TFBSTools::readJASPARMatrix("Data/ctcf_data/ctcf_motif.txt", matrixClass = "PFM")
  pwm <- TFBSTools::toPWM(pfm)

  chr_sels <- unique(tb_tads$chr)
  res_list <- vector("list", length = length(chr_sels))
  names(res_list) <- chr_sels
  for (chr_sel in chr_sels) {
    gr_ctcf <- gr_protein_gain("CTCF", chr_sel = chr_sel, chip_folder = chip_folder)
    res_list[[chr_sel]] <- tb_tads %>%
      dplyr::filter(chr %in% chr_sel) %>%
      calculate_CTCF_orientation_at_boundary(
        gr_ctcf = gr_ctcf, chr_sel = chr_sel, pwm = pwm
      )
  }

  rename_non_key <- function(df, keys, prefix) {
    dplyr::rename_with(df, ~paste0(prefix, .x), setdiff(names(df), keys))
  }

  tb_ctcf <- dplyr::bind_rows(res_list) %>%
    dplyr::filter(relScore >= 0.85) %>%
    dplyr::filter(!is.na(strand)) %>%
    dplyr::select(chr, tad_boundary, strand) %>%
    dplyr::group_by(chr, tad_boundary) %>%
    dplyr::summarise(strand = paste(sort(unique(strand)), collapse = ","), .groups = "drop") %>%
    dplyr::filter(!stringr::str_detect(strand, ","))
  ctcf_start <- rename_non_key(tb_ctcf, c("chr", "tad_boundary"), "start_")
  ctcf_end <- rename_non_key(tb_ctcf, c("chr", "tad_boundary"), "end_")

  tb_all <- tb_all %>%
    dplyr::mutate(start_p1 = start - 1) %>%
    dplyr::left_join(ctcf_start, by = c("chr", "start_p1" = "tad_boundary")) %>%
    dplyr::left_join(ctcf_end, by = c("chr", "end" = "tad_boundary")) %>%
    dplyr::select(-start_p1) %>%
    dplyr::mutate(
      orientation_different = dplyr::case_when(
        is.na(start_strand) | is.na(end_strand) ~ FALSE,
        stringr::str_detect(start_strand, "\\+") & stringr::str_detect(end_strand, "-") ~ TRUE,
        TRUE ~ FALSE
      )
    )

  res_part1 <- tb_all %>%
    dplyr::filter(!is.na(start_strand) & !is.na(end_strand)) %>%
    dplyr::select(meta.tool, meta.resol, meta.sub, chr, orientation_different) %>%
    dplyr::group_by(meta.tool, meta.resol, meta.sub) %>%
    dplyr::summarise(ctcf_prop = mean(orientation_different), .groups = "drop")
  res_part2 <- dplyr::select(tb_all, meta.tool, meta.resol, meta.sub) %>%
    dplyr::setdiff(dplyr::select(res_part1, meta.tool, meta.resol, meta.sub)) %>%
    dplyr::distinct()
  if (NROW(res_part2) > 0) {
    res <- res_part2 %>%
      dplyr::mutate(ctcf_prop = 0) %>%
      dplyr::bind_rows(res_part1)
  } else {
    res <- res_part1
  }

  res
}
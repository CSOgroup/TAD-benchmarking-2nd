source("utils/utils.R")

calculate_histone_bin_coverage <- function(gr_bin, gr_histone, binsize) {
  df_intersect <- GenomicRanges::findOverlaps(
    gr_bin, gr_histone,
    minoverlap = 2
  ) %>%
    tibble::as_tibble()
  part1.a <- gr_bin[df_intersect$queryHits]
  part1.b <- gr_histone[df_intersect$subjectHits]
  part1.width <- GenomicRanges::pintersect(part1.a, part1.b) %>%
    GenomicRanges::width() %>%
    magrittr::subtract(1)
  part1.a %>%
    tibble::as_tibble() %>%
    dplyr::select(chr = seqnames, start, end) %>%
    dplyr::mutate(fc_bp = part1.b$vlaue * part1.width) %>%
    dplyr::group_by(chr, start, end) %>%
    dplyr::summarise(fc_bp = sum(fc_bp), .groups = "drop") %>%
    dplyr::mutate(fc_bp = fc_bp / binsize) %>%
    dplyr::filter(fc_bp > 0)
}


prepare_chromosome_histone_data <- function(binsize, gr_h27 = NULL, gr_h36 = NULL, chr_sel = "chr6") {
  `%>%` <- magrittr::`%>%`
  chr_len <- chrsize_hg19()[[chr_sel]]
  if (chr_len %% binsize > 0) {
    gr_binchr <- tibble::tibble(
      chr = chr_sel,
      start = seq(1, chr_len, binsize),
      end = c(seq(binsize, chr_len, binsize), chr_len)
    ) %>%
      GenomicRanges::GRanges()
  } else {
    gr_binchr <- tibble::tibble(
      chr = chr_sel,
      start = seq(1, chr_len - binsize, binsize),
      end = seq(binsize, chr_len, binsize)
    ) %>%
      GenomicRanges::GRanges()
  }
  if (is.null(gr_h27)) {
    gr_h27 <- gr_protein_gain("H3K27me3", chr_sel = chr_sel)
  }
  if (is.null(gr_h36)) {
    gr_h36 <- gr_protein_gain("H3K36me3", chr_sel = chr_sel)
  }
  h27_binned <- calculate_histone_bin_coverage(gr_binchr, gr_h27, binsize)
  h36_binned <- calculate_histone_bin_coverage(gr_binchr, gr_h36, binsize)
  colnames(h27_binned) <- c("chr", "start", "end", "reads_h27")
  colnames(h36_binned) <- c("chr", "start", "end", "reads_h36")
  h_binned <- dplyr::inner_join(
    h27_binned, h36_binned, by = c("chr", "start", "end")
  ) %>%
    dplyr::filter(reads_h27 > 0, reads_h36 > 0) %>%
    dplyr::mutate(lr = log10(reads_h27 / reads_h36))
  h_binned %>%
    GenomicRanges::GRanges()
}

generate_shuffled_histone_data <- function(gr_h_binned, nshuf = 10) {
  nrow_h_binned <- length(gr_h_binned)
  res <- vector("list", nshuf)
  for (i in seq(nshuf)) {
    shuf_pos <- sample(seq(nrow_h_binned), size = nrow_h_binned)
    gr_h_binned$lr <- gr_h_binned$lr[shuf_pos]
    res[[i]] <- gr_h_binned
  }
  res
}


calculate_log_ratio_histone <- function(gr_tads, gr_h_binned) {
  df_intersect <- GenomicRanges::findOverlaps(gr_tads, gr_h_binned) %>%
    tibble::as_tibble()
  if (NROW(df_intersect) == 0) {
    return(rep(0, length(gr_tads)))
  } else {
    part1.a <- gr_tads[df_intersect$queryHits]
    part1.b <- gr_h_binned[df_intersect$subjectHits]
    part1.width <- GenomicRanges::pintersect(part1.a, part1.b) %>%
      GenomicRanges::width() %>%
      magrittr::subtract(1)
    lr_res <- part1.a %>%
      tibble::as_tibble() %>%
      dplyr::mutate(lr_part = part1.width * part1.b$lr) %>%
      dplyr::group_by(seqnames, start, end) %>%
      dplyr::summarise(lr_sum = sum(lr_part), .groups = "drop") %>%
      dplyr::mutate(lr_mean = lr_sum / (start - end))
    res <- c(lr_res$lr_mean, rep(0, length(gr_tads) - NROW(lr_res)))
  }
  res
}


analyze_histone_modification_single_tool_genome <- function(tb_tads, gr_h27 = NULL, gr_h36 = NULL,
                                                            share = 0.1, fdr_thresh = 0.1,
                                                            nshuf = 10, chr_sel = "chr6") {
  gr_tads <- tb_tads %>%
    GenomicRanges::makeGRangesFromDataFrame()
  binsize <- round(mean(GenomicRanges::width(gr_tads) - 1) * share)
  gr_h_binned <- prepare_chromosome_histone_data(binsize, gr_h27 = gr_h27, gr_h36 = gr_h36, chr_sel = chr_sel)
  gr_h_binned_list <- generate_shuffled_histone_data(gr_h_binned, nshuf = nshuf)
  avg_lr <- calculate_log_ratio_histone(gr_tads, gr_h_binned)
  distr_shuf_lr <- lapply(
    gr_h_binned_list, calculate_log_ratio_histone, gr_tads = gr_tads
  ) %>%
    unlist()
  median_shuf_lr <- median(distr_shuf_lr)
  distr_shuf_lr_len <- length(distr_shuf_lr)
  pval <- vector(length = length(avg_lr))
  for (i in seq_along(avg_lr)) {
    query <- avg_lr[i]
    if (query >= median_shuf_lr) {
      pval[i] <- sum(distr_shuf_lr >= query) / distr_shuf_lr_len
    } else {
      pval[i] <- sum(distr_shuf_lr <= query) / distr_shuf_lr_len
    }
  }
  qval_BH <- p.adjust(pval, method = "BH")
  # qval_Bonferroni <- p.adjust(pval, method = "bonferroni")
  shareTadsSignif_lr_BH_sum <- sum(qval_BH < fdr_thresh)
  tads_n <- nrow(tb_tads)

  tb_tads %>%
    dplyr::select(meta.tool, meta.resol, meta.sub, chr) %>%
    dplyr::distinct() %>%
    dplyr::mutate(shareTadsSignif_lr_BH_sum = shareTadsSignif_lr_BH_sum, tads_n = tads_n)
}


analyze_histone_modification_all_tools_genome <- function(tb_tads, gr_h27 = NULL, gr_h36 = NULL,
                                                          share = 0.1, fdr_thresh = 0.1,
                                                          nshuf = 10, chr_sel = "chr6") {
  input_list <- tb_tads %>%
    dplyr::mutate(grp = paste(meta.tool, meta.resol, meta.sub, chr, sep = "||")) %>%
    split(.$grp)

  res <- input_list %>%
    lapply(
      analyze_histone_modification_single_tool_genome, gr_h27 = gr_h27, gr_h36 = gr_h36,
      share = share, fdr_thresh = fdr_thresh, nshuf = nshuf,
      chr_sel = chr_sel) %>%
    dplyr::bind_rows()

  res
}

analyze_histone_modification_all_tools <- function(tb_tads, chip_folder = "Data") {
  # chr_sels <- paste0("chr", c(1:22, "X"))
  chr_sels <- unique(tb_tads$chr)
  res_list <- vector("list", length = length(chr_sels))
  names(res_list) <- chr_sels
  for (chr_sel in chr_sels) {
    gr_h27 <- gr_protein_gain("H3K27me3", chr_sel = chr_sel, chip_folder = chip_folder)
    gr_h36 <- gr_protein_gain("H3K36me3", chr_sel = chr_sel, chip_folder = chip_folder)
    res_list[[chr_sel]] <- tb_tads %>%
      dplyr::filter(chr %in% chr_sel) %>%
      analyze_histone_modification_all_tools_genome(gr_h27 = gr_h27, gr_h36 = gr_h36,
                                                    share = 0.1, fdr_thresh = 0.1,
                                                    nshuf = 10, chr_sel = chr_sel)
  }

  dplyr::bind_rows(res_list) %>%
    dplyr::group_by(meta.tool, meta.resol, meta.sub) %>%
    dplyr::summarise(shareTadsSignif_lr_BH_sum = sum(shareTadsSignif_lr_BH_sum), tads_n = sum(tads_n), .groups = "drop") %>%
    dplyr::mutate(shareTadsSignif_lr_BH = shareTadsSignif_lr_BH_sum / tads_n) %>%
    dplyr::select(-shareTadsSignif_lr_BH_sum, -tads_n)
}

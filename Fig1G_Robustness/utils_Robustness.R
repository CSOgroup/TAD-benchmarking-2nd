source("utils/utils.R")


moc_overlap_widthfirst <- function(gr) {
  gr <- sort(gr)
  gr$width_index <- GenomicRanges::width(gr)

  gr_sel <- function(gr) {
    res <- gr[gr$width_index == max(gr$width_index)]
    if (length(res) != 1) {
      res <- res[1]
    }
    list(
      res = res,
      gr_out = gr[!IRanges::overlapsAny(gr, res)]
    )
  }

  res <- vector("list", length(gr))
  first_index <- IRanges::countOverlaps(gr, gr)
  res[[1]] <- list(
    res = gr[first_index == 1],
    gr_out = gr[first_index > 1]
  )
  # res[[1]] <- gr_sel(gr)
  for (i in seq(2, length(res))) {
    if (length(res[[i - 1]]$gr_out) == 0) {
      break
    }
    res[[i]] <- gr_sel(res[[i - 1]]$gr_out)
  }
  res_df <- res %>%
    lapply(function(x) { tibble::as_tibble(x$res) }) %>%
    dplyr::bind_rows() %>%
    dplyr::rename(chr = seqnames) %>%
    dplyr::select(-c(width, strand, width_index))
  # dplyr::select(chr = seqnames, start, end)
  res_df
}

moc_fill_part <- function(fillDT, chrSize, chr_sel = "chr6") {
  colnames(fillDT) <- c("chr", "start", "end")
  gr_whole <- GenomicRanges::GRanges(
    tibble::tibble(chr = chr_sel, start = 1, end = chrSize)
  )
  gr_tads <- GenomicRanges::GRanges(fillDT)
  if (any(IRanges::countOverlaps(gr_tads, gr_tads) > 1)) {
    gr_tads <- moc_overlap_widthfirst(gr_tads) %>%
      GenomicRanges::GRanges() %>%
      sort()
  }
  gr_inter <- GenomicRanges::setdiff(gr_whole, gr_tads)
  c(gr_tads, gr_inter) %>%
    sort()
}


moc_fromranges <- function(set1DT, set2DT) {
  if (NROW(set1DT) == 0 && NROW(set2DT) > 0) {
    MoC_score <- 0
  } else if (NROW(set1DT) > 0 && NROW(set2DT) == 0) {
    MoC_score <- 0
  } else if ((NROW(set1DT) == NROW(set2DT)) &&
    all(GenomicRanges::start(set1DT) == GenomicRanges::start(set2DT)) &&
    all(GenomicRanges::end(set1DT) == GenomicRanges::end(set2DT))) {
    MoC_score <- 1
  } else {
    df_intersect <- GenomicRanges::findOverlaps(set1DT, set2DT) %>%
      tibble::as_tibble()
    if (NROW(df_intersect) == 0) {
      MoC_score <- 0
    } else {
      set1DT.long <- set1DT[df_intersect$queryHits]
      set2DT.long <- set2DT[df_intersect$subjectHits]
      intersect.width <- GenomicRanges::pintersect(set1DT.long, set2DT.long) %>%
        GenomicRanges::width() %>%
        as.double()
      set1DT.long.width <- set1DT.long %>%
        GenomicRanges::width() %>%
        as.double()
      set2DT.long.width <- set2DT.long %>%
        GenomicRanges::width() %>%
        as.double()
      MoC_score <- (
        sum(intersect.width^2 / (set1DT.long.width * set2DT.long.width)) - 1
      ) / (sqrt(NROW(set1DT) * NROW(set2DT)) - 1)
    }
  }
  MoC_score
}


moc_singletool <- function(tb_tads_onetool, chrSize, chr_sel = "chr6") {
  subs <- tb_tads_onetool$meta.sub %>% setdiff("sub1")
  tools_new <- unique(tb_tads_onetool$meta.tool)
  set2DT <- tb_tads_onetool %>%
    dplyr::filter(meta.sub == "sub1") %>%
    dplyr::select(chr, start, end) %>%
    moc_fill_part(chrSize, chr_sel = chr_sel)

  sub_dfs <- tb_tads_onetool %>%
    dplyr::filter(meta.sub %in% subs)
  sub_ranges <- split(sub_dfs, sub_dfs$meta.sub) %>%
    lapply(dplyr::select, chr, start, end) %>%
    lapply(moc_fill_part, chrSize, chr_sel = chr_sel) %>%
    lapply(moc_fromranges, set2DT) %>%
    unlist()
  sub1 <- 1
  names(sub1) <- "sub1"
  sub_all <- c(sub_ranges, sub1)
  tibble::tibble(
    meta.tool = tools_new,
    meta.sub = names(sub_all),
    value = sub_all
  )
}


moc_multitool <- function(tb_tads, chr_sel = "chr6") {
  chrSize <- chrsize_hg19()[[chr_sel]]
  input_list <- tb_tads %>%
    dplyr::arrange(meta.tool, meta.resol, meta.sub, start) %>%
    split(~meta.tool) %>%
    purrr::discard(~NROW(.) == 0)
  res <- input_list %>%
    lapply(moc_singletool, chrSize = chrSize, chr_sel = chr_sel) %>%
    dplyr::bind_rows()
  res
}


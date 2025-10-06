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

moc_fill_part <- function(fillDT) {
  colnames(fillDT) <- c("chr", "start", "end")
  gr_whole <- GenomicRanges::GRanges(
    tibble::tibble(
      chr = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX"),
      start = 1,
      end = c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560)
    )
  )
  gr_tads <- sort(GenomicRanges::GRanges(fillDT))
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
      ) / (sqrt(as.double(NROW(set1DT)) * as.double(NROW(set2DT))) - 1)
    }
  }
  MoC_score
}


moc_singletool <- function(tb_tads_onetool) {
  subs <- tb_tads_onetool$meta.sub %>% setdiff("sub1")
  tools_new <- unique(tb_tads_onetool$meta.tool)
  meta.resol <- unique(tb_tads_onetool$meta.resol)
  if (length(subs) == 0) {
    return(tibble::tibble(
      meta.tool = tools_new,
      meta.resol = meta.resol,
      meta.sub = c("sub0.0001", "sub0.0002", "sub0.0005", "sub0.001", "sub0.002", "sub0.005", "sub0.01", "sub0.02", "sub0.05", "sub0.1", "sub0.2", "sub0.5", "sub1"),
      value = 0
    ))
  }

  sub1_df <- tb_tads_onetool %>%
    dplyr::filter(meta.sub == "sub1")
  if (NROW(sub1_df) == 0) {
    return(tibble::tibble(
      meta.tool = tools_new,
      meta.resol = meta.resol,
      meta.sub = c("sub0.0001", "sub0.0002", "sub0.0005", "sub0.001", "sub0.002", "sub0.005", "sub0.01", "sub0.02", "sub0.05", "sub0.1", "sub0.2", "sub0.5", "sub1"),
      value = 0
    ))
  }
  set2DT <- sub1_df %>%
    dplyr::select(chr, start, end) %>%
    moc_fill_part()

  sub_dfs <- tb_tads_onetool %>%
    dplyr::filter(meta.sub %in% subs)
  sub_ranges <- split(sub_dfs, sub_dfs$meta.sub) %>%
    lapply(function(x) {
      dplyr::select(x, chr, start, end) %>%
        moc_fill_part() %>%
        moc_fromranges(set2DT)
    }) %>%
    unlist()
  sub1 <- 1
  names(sub1) <- "sub1"
  sub_all <- c(sub_ranges, sub1)
  tibble::tibble(
    meta.tool = tools_new,
    meta.resol = meta.resol,
    meta.sub = names(sub_all),
    value = sub_all
  )
}


moc_multitool <- function(tb_tads) {
  input_list <- tb_tads %>%
    dplyr::arrange(meta.tool, meta.resol, meta.sub, start) %>%
    split(~meta.tool) %>%
    purrr::discard(~NROW(.) == 0)
  res <- input_list %>%
    lapply(moc_singletool) %>%
    dplyr::bind_rows()
  res
}


source("utils/utils.R")

calculate_insulation_scores <- function(adj, chr_sel = "chr6", bin = 50) {
  chrSize <- chrsize_hg19()[[chr_sel]]
  row_max <- ceiling(chrSize / 10000)
  padding_n <- row_max - nrow(adj)
  if (padding_n != 0) {
    adj_padding <- cbind(adj, matrix(0, nrow = nrow(adj), ncol = padding_n)) %>%
      rbind(matrix(0, nrow = padding_n, ncol = padding_n + nrow(adj)))
  } else {
    adj_padding <- adj
  }


  calculate_boundary_score <- function(boundary, mx, bin) {
    mx[(boundary - bin):(boundary - 1), ((boundary + 1):(boundary + bin))] %>%
      mean()
  }


  n <- nrow(adj_padding)
  mx_boundaries <- seq_len(n)
  mx_boundaries2 <- mx_boundaries[
    mx_boundaries - bin >= 1 & mx_boundaries + bin <= n
  ]
  boundary_is <- mx_boundaries2 %>%
    lapply(calculate_boundary_score, mx = adj_padding, bin = bin) %>%
    unlist()
  tibble::tibble(boundary = mx_boundaries2, value = boundary_is)
}


calculate_boundary_insulation_ratio <- function(boundaries, df_ISindex, gap_len = 3, bg_len = 30) {
  boundaries_low_cut <- min(df_ISindex$boundary) + gap_len + bg_len
  boundaries_high_cut <- max(df_ISindex$boundary) - gap_len - bg_len
  boundaries2 <- boundaries[boundaries >= boundaries_low_cut &
                              boundaries <= boundaries_high_cut]

  boundaries2_bg <- boundaries2 %>%
    lapply(function(x) { setdiff(c(
      seq(x - gap_len - bg_len, x - gap_len - 1),
      seq(x + gap_len + 1, x + gap_len + bg_len)
    ), boundaries2) })
  boundary_is <- df_ISindex$value
  names(boundary_is) <- df_ISindex$boundary
  pos <- boundary_is[as.character(boundaries2)]
  neg <- boundaries2_bg %>%
    lapply(function(x) { mean(boundary_is[as.character(x)]) }) %>%
    unlist()
  mean_res <- pos / neg
  res <- mean(mean_res[!is.infinite(mean_res)], na.rm = TRUE)

  res
}


calculate_insulation_scores_for_multiple_tools <- function(df_diffresol, gap_len = 3, bg_len = 30, df_ISindex = NULL) {
  resols <- unique(df_diffresol$meta.resol)
  res <- vector("list", length = length(resols))
  for (i in seq_along(res)) {
    resol <- resols[i]
    df_diffresol_list <- df_diffresol %>%
      dplyr::filter(meta.resol %in% resol) %>%
      split(.$meta.tool) %>%
      lapply(function(df) {
        # unique(round(c((df$start - 1) / 10000 + 1, df$end / 10000)))
        df %>%
          dplyr::mutate(startbin = round((start - 1) / 10000 + 1),
                        endbin = round(end / 10000))
      })

    res[[i]] <- df_diffresol_list %>%
      lapply(function(df) (unique(c(df$startbin, df$endbin)))) %>%
      lapply(calculate_boundary_insulation_ratio, df_ISindex, gap_len = gap_len, bg_len = bg_len) %>%
      list2tb() %>%
      dplyr::arrange(lt_value) %>%
      dplyr::mutate(meta.resol = resol) %>%
      dplyr::select(meta.tool = lt_name, meta.resol, is = lt_value)
  }
  dplyr::bind_rows(res) %>%
    dplyr::arrange(is)
}

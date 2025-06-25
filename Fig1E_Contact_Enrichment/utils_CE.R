source("utils/utils.R")

calculate_contact_enrichment_matrix <- function(adj, chr_sel = "chr6") {
  chrSize <- chrsize_hg19()[[chr_sel]]
  row_max <- ceiling(chrSize / 10000)
  padding_n <- row_max - nrow(adj)
  if (padding_n != 0) {
    adj_padding <- cbind(adj, matrix(0, nrow = nrow(adj), ncol = padding_n)) %>%
      rbind(matrix(0, nrow = padding_n, ncol = padding_n + nrow(adj)))
  } else {
    adj_padding <- adj
  }

  n <- nrow(adj_padding)
  mx_index <- abs(
    matrix(rep(seq(n), n), n) - matrix(rep(seq(n), n), n, byrow = TRUE)
  )
  replace_num <- adj_padding %>%
    split(mx_index) %>%
    lapply(mean)
  mx_expcted <- plyr::mapvalues(
    mx_index,
    as.double(names(replace_num)),
    unlist(replace_num)
  )
  ((adj_padding + 1) / (mx_expcted + 1))
}


calculate_tad_contact_enrichment <- function(tb_tads, enrichment_matrix = NULL) {

  expandplus <- function(df) {
    a <- seq(df$binstart, df$binend)
    expand.grid(X1 = a, X2 = a) %>%
      dplyr::mutate(start = df$start, end = df$end)
  }

  tb_input2 <- tb_tads %>%
    dplyr::mutate(
      binstart = floor((start - 1) / 10000) + 1,
      binend = floor(end / 10000),
      index = seq(NROW(.))
    )
  modularity_con <- tb_input2 %>%
    split(.$index) %>%
    lapply(function(df) {
      a <- seq(df$binstart, df$binend)
      expand.grid(X1 = a, X2 = a)
    }) %>%
    dplyr::bind_rows()

  modularity_con <- dplyr::distinct(modularity_con)

  enrichment_matrix_n <- NROW(enrichment_matrix)
  modularity_con <- modularity_con %>%
    dplyr::filter(X1 <= enrichment_matrix_n, X2 <= enrichment_matrix_n)
  res <- mean(enrichment_matrix[as.matrix(modularity_con[, 1:2])])

  res
}

calculate_contact_enrichment_for_multiple_tools <- function(df_diffresol, enrichment_matrix = NULL) {
  df_diffresol_list <- df_diffresol %>%
    split(.$meta.tool)
  if (length(df_diffresol_list) == 1) {
    enrichment_scores <- calculate_tad_contact_enrichment(
      df_diffresol_list[[1]],
      enrichment_matrix = enrichment_matrix
    )
  } else {
    enrichment_scores <- lapply(
      df_diffresol_list, calculate_tad_contact_enrichment,
      enrichment_matrix = enrichment_matrix
    )
  }

  tibble::tibble(
    meta.tool = names(df_diffresol_list),
    contact_enrichment = unlist(enrichment_scores)
  ) %>%
    dplyr::left_join(
      dplyr::distinct(dplyr::select(df_diffresol, meta.tool, meta.resol)),
      by = "meta.tool"
    ) %>%
    dplyr::select(meta.tool, meta.resol, contact_enrichment)
}

`%>%` <- magrittr::`%>%`

gr_protein_gain <- function(protein, chr_sel = "chr6", chip_folder = "Data") {
  if (protein %in% c("CTCF", "RAD21", "SMC3")) {
    file_path <- file.path(chip_folder, sprintf("ProteinPeaks/%s_peaks.bed", protein))

    res <- file_path %>%
      readr::read_tsv(col_names = FALSE, show_col_types = FALSE) %>%
      dplyr::filter(X1 == chr_sel) %>%
      dplyr::select(chr = X1, start = X2, end = X3) %>%
      GenomicRanges::GRanges()
  } else if (protein %in% c("H3K27me3", "H3K36me3")) {
    filename <- ifelse(
      protein == "H3K27me3",
      sprintf("ENCFF594HSG_%s.bedGraph.gz", chr_sel),  # H3K27me3 - repression
      sprintf("ENCFF662QFK_%s.bedGraph.gz", chr_sel)   # H3K36me3 - active
    )

    file_path <- file.path(chip_folder, "HistoneSignals", filename)

    res <- file_path %>%
      readr::read_tsv(
        col_names = c("chr", "start", "end", "vlaue"),
        show_col_types = FALSE, comment = "#"
      ) %>%
      GenomicRanges::GRanges()
  }
  return(res)
}

chrsize_hg19 <- function() {
  res <- list(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560, 59373566)
  names(res) <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")
  res
}

bin_trim <- function(gr, small_bin) {

  GenomicRanges::width(gr) <- floor(
    GenomicRanges::width(gr) / small_bin
  ) * small_bin
  GenomicRanges::tile(gr[GenomicRanges::width(gr) > 0], width = small_bin)
}


list2tb <- function(lt, con = c("s", "df")) {

  con <- con[1]
  if (con == "s") {
    list_names <- names(lt)
    res <- vector("list", length = length(lt))
    for (i in seq(lt)) {
      res[[i]] <- tibble::tibble(
        lt_name = list_names[i],
        lt_value = lt[[i]]
      )
    }
    dplyr::bind_rows(res)
  } else if (con == "df") {
    lt_name <- rep(names(lt), lapply(lt, NROW))
    dplyr::bind_rows(lt) %>%
      dplyr::mutate(lt_name = lt_name) %>%
      dplyr::select(lt_name, dplyr::everything())
  }
}

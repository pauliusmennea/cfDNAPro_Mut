#' plot copy number profile
#' @import ggplot2
#' @import tibble
#' @import dplyr
#' @import TxDb.Hsapiens.UCSC.hg38.knownGene
#' @import TxDb.Hsapiens.UCSC.hg19.knownGene
#' @import Homo.sapiens
#' @importFrom rlang has_name
#' @importFrom GenomicFeatures genes
#' @importFrom ggrepel geom_text_repel
#' @importFrom purrr map
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom plyranges join_overlap_inner expand_ranges
#'
#' @param x QDNAseqCopyNumbers object
#' @param mut_gr optional GRanges fragment/mutation object to overlap with CN
#' @param gene_to_highlight A named list.
#' @param genome string. hg19 or hg38
#' @param ylim vector. default is c(-2, 2).
#' @param chromosome vector. default is c(seq(1, 22, 1), "X").
#' @param point_color named vector
#' @param x_title string.
#' @param y_title string.
#' @param point_size numerical. default is 0.3.
#' @param point_alpha numerical. default is 0.9.
#' @param chr_edge_color string. default is "black".
#' @param chr_edge_line_size numerical. default is 0.2.
#' @param chr_edge_alpha numerical. default is 0.8.
#' @param chr_edge_type string. default is "dotted".
#' @param segment_color string. default is "red".
#' @param segment_alpha numerical. default is 1.
#' @param segment_line_end string. default is "round".
#' @param segment_line_size numerical. default is 0.75.
#' @param legend_position string. default is "none", which mean no legends.
#' @param x_axis_expand numerical vector. default is c(0.1, 0.1).
#' @param y_axis_expand numerical vector. default is c(0, 0).
#' @param text_size numerical. relevant when mut_gr is provided. default is 1.6.
#' @param ... Other params to geom_txt_repel()
#'
#' @return This function returns ggplot2 object.
#' @export
#' @author Haichao Wang
#'
#' @examples
#' \dontrun{
#'
#' p <- plotCNV(bamfile = "/path/to/bamfile.bam")
#' }
plotCNV <- function(x, mut_gr,
                    gene_to_highlight = list("ENTREZID" = NULL,
                                             "ENSEMBL" = NULL,
                                             "SYMBOL" = c("TTC34")),
                    genome = "hg19",
                    ylim = c(-2, 2),
                    chromosome = c(seq(1, 22, 1), "X"),
                    point_color = c("Loss" = "royalblue",
                                    "Deletion" = "darkblue",
                                    "Neutral" = "darkgrey",
                                    "Amplification" = "darkorange",
                                    "Gain" = "orange3"),
                    x_title = "Chromosome",
                    y_title = "Log2 Ratio",
                    point_size = 0.3,
                    point_alpha = 0.9,
                    chr_edge_color = "black",
                    chr_edge_line_size = 0.2,
                    chr_edge_alpha = 0.8,
                    chr_edge_type = "dotted",
                    segment_color = "red",
                    segment_alpha = 1,
                    segment_line_end = "round",
                    segment_line_size = 0.75,
                    legend_position = "none",
                    x_axis_expand = c(0.1, 0.1),
                    y_axis_expand = c(0, 0),
                    text_size = 1.6,
                    ...) {

  names(point_color) <- stringr::str_to_title(names(point_color))

  if (!rlang::has_name(point_color, "Loss")) {
    point_color["Loss"] <- "royalblue"
  }
  if (!rlang::has_name(point_color, "Deletion")) {
    point_color["Deletion"] <- "darkblue"
  }
  if (!rlang::has_name(point_color, "Neutral")) {
    point_color["Neutral"] <- "darkgrey"
  }
  if (!rlang::has_name(point_color, "Amplification")) {
    point_color["Amplification"] <- "darkorange"
  }
  if (!rlang::has_name(point_color, "Gain")) {
    point_color["Gain"] <- "orange3"
  }

  # extract QDNAseqCopyNumbers data
  if (is(x, "QDNAseqCopyNumbers")) {

    sample_name <- x@phenoData@data$name

    chr_levels <- as.character(c(seq(1, 22, 1), "X", "Y"))

    cnv_tibble <- tibble::tibble(
      chr = x@featureData@data$chromosome,
      start = x@featureData@data$start,
      end = x@featureData@data$end,
      use = x@featureData@data$use,
      copynumber = x@assayData$copynumber[, sample_name],
      seg = x@assayData$segmented[, sample_name],
      call = x@assayData$calls[, sample_name]
    )

    cnv_tibble$chr <- factor(cnv_tibble$chr, levels = chr_levels)

    cnv_tibble2 <- cnv_tibble %>%
      dplyr::filter(.data$chr %in% chromosome) %>%
      dplyr::filter(.data$use == TRUE) %>%
      dplyr::filter(is.finite(.data$copynumber)) %>%
      dplyr::mutate(call = dplyr::case_when(
        .data$call == -2 ~ "Deletion",
        .data$call == -1 ~ "Loss",
        .data$call == 0 ~ "Neutral",
        .data$call == 1 ~ "Gain",
        .data$call == 2 ~ "Amplification"
      )) %>%
      dplyr::mutate(x_index = dplyr::row_number()) %>%
      dplyr::mutate(strand = "*") %>%
      dplyr::mutate(logratio = log2(.data$copynumber)) %>%
      dplyr::mutate(logseg = log2(.data$seg)) %>%
      dplyr::mutate(seg_id = rleid_cfdnapro(.data$logseg))

    cnv_tibble2$call <- factor(cnv_tibble2$call,
                               levels = c("Amplification", "Gain", "Neutral",
                                          "Loss", "Deletion"))

    cnv_tibble3 <- cnv_tibble2 %>%
      dplyr::mutate(seqnames = paste("chr", chr, sep = ""))

    cnv_tibble3_gr <- makeGRangesFromDataFrame(
      df = cnv_tibble3,
      seqnames.field = "seqnames",
      keep.extra.columns = TRUE
    )

  } else {
    stop("Only QDNAseqCopyNumbers object is accepted as input for the moment. ")
  }

  # obtain gene position
  if (genome == "hg38") {
    message("Switching to UCSC hg38 reference genome.")
    TxDb(Homo.sapiens) <- TxDb.Hsapiens.UCSC.hg38.knownGene
  }

  hs <- Homo.sapiens
  all_genes <- suppressMessages(genes(hs,
                                      columns = c("ENTREZID",
                                                  "ENSEMBL",
                                                  "SYMBOL")))

  expanded_all_genes <- expand_ranges(all_genes,
                                      "ENSEMBL",
                                      .keep_empty = TRUE) %>%
    expand_ranges("SYMBOL", .keep_empty = TRUE) %>%
    expand_ranges("ENTREZID", .keep_empty = TRUE)

  # filter by targeted regions/genes
  expanded_all_genes_tb <- tibble::as_tibble(expanded_all_genes)
  ensembl <- gene_to_highlight[["ENSEMBL"]]
  symbol <- gene_to_highlight[["SYMBOL"]]
  entrezid <- gene_to_highlight[["ENTREZID"]]

  filter_fun <- function(x){
    col <- names(gene_to_highlight)[[x]]

    if(is.null(gene_to_highlight[[x]])) {

      return(NULL)
    } else {
      ans <- plyranges::filter(expanded_all_genes_tb,
                               .data[[col]] %in% gene_to_highlight[[x]])
      message("Found these entries: \n")
      print(ans)

      items_not_found <- setdiff(gene_to_highlight[[x]], ans[[col]])

      if (length(items_not_found) != 0) {
        message("These entries are not found: \n")
        message(items_not_found)
      }

      return(ans)
    }
  }

  filter_ans <- map(seq_len(length(gene_to_highlight)), .f = filter_fun) %>%
    dplyr::bind_rows() %>%
    dplyr::select(-width)

  filter_ans_gr <- makeGRangesFromDataFrame(df = filter_ans,
    keep.extra.columns = TRUE
  )

  # overlap bin and gene
  olap <- plyranges::join_overlap_inner(cnv_tibble3_gr, filter_ans_gr)

  # if olap also overlaps GRanges mutations, then add a mut column
  if (!is.null(mut_gr)) {

    olap_mut <- plyranges::join_overlap_inner(olap, mut_gr)

    olap_mut_df <- as.data.frame(olap_mut)
    olap_mut_df <- olap_mut_df %>% dplyr::select(c(1:17, 21))

    olap_df <- olap_mut_df %>%
      group_by_at(vars(1:(ncol(olap_mut_df) - 1))) %>%
      summarize(MUT_count = sum(grepl("MUT:concordant|MUT:single_read",
                locus_status)),
                ALL_count = sum(grepl("[A-Za-z]", locus_status))) %>%
      ungroup() %>%
      as.data.frame()

  } else {

    olap_df <- as.data.frame(olap)

  }

  # create chromosome boundaries
  chr_edge <- cnv_tibble2 %>%
    dplyr::group_by(chr) %>%
    dplyr::filter(row_number() == n()) %>%
    dplyr::pull(x_index)

  chr_edge <- c(1, chr_edge)

  x_breaks <- cnv_tibble2 %>%
    dplyr::group_by(chr) %>%
    dplyr::filter(row_number() == ceiling(n()/2)) %>%
    dplyr::pull(x_index)

  x_labels <- as.character(unique(cnv_tibble2$chr))

  # ggplot2 plot
  p <- ggplot(data = cnv_tibble2) +
    # copy number points
    geom_point(aes(x = .data$x_index, y = .data$logratio, color = .data$call),
               size = point_size,
               alpha = point_alpha) +
    scale_color_manual(values = point_color) +
    # segmentation line
    geom_line(aes(x = .data$x_index, y = .data$logseg, group = .data$seg_id),
              colour = segment_color,
              alpha = segment_alpha,
              size = segment_line_size,
              lineend = segment_line_end) +
    # chr edges
    geom_vline(xintercept = chr_edge,
               linetype = chr_edge_type,
               color = chr_edge_color,
               size = chr_edge_line_size) +
    # chr names
    scale_x_continuous(
      breaks = x_breaks,
      labels = x_labels,
      expand = c(0.01, 0.01)) +
    scale_y_continuous(limits = ylim,
                       expand = y_axis_expand) +
    labs(x = x_title, y = y_title) +
    theme_cnv_plot()  +
    theme(legend.position = legend_position)

  # check if 'MUT_count' column exists in olap_df
  # and add MUT/REF fragment level annotations to the plot
  if ("MUT_count" %in% names(olap_df)) {
    p <- p +
      geom_point(data = olap_df,
                mapping = aes(x = .data$x_index, y = .data$logratio),
                size = point_size + 1,
                color = "black",
                shape = 1
      ) +
      geom_text_repel(data = olap_df,
                      aes(x = .data$x_index,
                          y = .data$logratio,
                          label = paste(.data$SYMBOL, " Gene: ",
                                        MUT_count, "/",
                                        ALL_count,
                          " Fragments Carrying Candidate Mutations",
                          sep = "")),
                      box.padding = 0.5,
                      min.segment.length = 0,
                      segment.linetype = 1,
                      segment.size = segment_line_size / 2,
                      segment.color = "black",
                      nudge_y = 0.3,
                      max.overlaps = Inf,
                      size = text_size
      )
    print("CN plot with integrated mutational information.")
    print("The GRanges object contains the mutation data.")
    print("Selected genes are annotated with fragment counts.")
    print("Only pre-specified SNVs are used to generate these counts.")
    print("The SNV loci can be indicated in readBam().")
    print("All SNVs within a gene are reflected in the fragment ratio.")
  } else {
    print("Default CN plot.")
  }
  return(p)
}

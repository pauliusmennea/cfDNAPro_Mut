
#' call fragment length
#' @import ggplot2
#' @import dplyr
#' @import plyranges
#'
#' @param fragment_obj
#' @param fragment_type String; "default", "mut_vs_ref", "mut_vs_outer", "mut_vs_ref_norm", "mut_vs_outer_norm"
#' @param isize_min
#' @param isize_max
#' @param ...
#'
#' @return a tibble
#' @export
#'
#' @examples
callLength <- function(fragment_obj,
                       fragment_type = "default",
                       isize_min = 1L,
                       isize_max = 1000L,
                       ...) {
  if (fragment_type == "default") {
    frag <- fragment_obj
    # calculating insert sizes
    message("Calculating insert sizes...")
    frag$insert_size <- BiocGenerics::width(frag)

    # size analysis
    frag <- plyranges::filter(frag,
                              insert_size >= isize_min & insert_size <= isize_max)
    isize <- frag$insert_size

    isize_tibble <- tibble("insert_size" = isize, "count" = 1) %>%
      dplyr::filter(!is.na(insert_size))

    result <- isize_tibble %>%
      dplyr::group_by(.data$insert_size) %>%
      dplyr::summarise("All_Reads.fr_count" = sum(count))

    # quality control results
    # Create a vector of elements
    isize_ref <- seq.int(isize_min, isize_max, by = 1L) %>%
      as_tibble()

    colnames(isize_ref) <- c("insert_size")

    # report abnormal isizes

    missing_isize <- dplyr::anti_join(isize_ref, result, by = "insert_size")

    # handle missing isize(s)

    if (nrow(missing_isize) != 0) {
      message("Missing isize detected: ")
      print(dplyr::pull(missing_isize, insert_size))

      result <- dplyr::right_join(result, isize_ref, by = "insert_size") %>%
        tidyr::replace_na(replace = list(All_Reads.fr_count = 0)) %>%
        dplyr::arrange(insert_size) %>%
        dplyr::mutate(prop = All_Reads.fr_count / sum(All_Reads.fr_count))

      message("Missing isize(s) added back to the final result with count of 0!")

      message("Job completed successfully. ")
    }

  } else if (grepl("mut", fragment_type, ignore.case = TRUE)) {
    # Convert fragment_obj to data frame
    # and filter out discordant and non-MUT base fragments
    size_table_inner <- as.data.frame(fragment_obj) %>%
      filter(!grepl("discordant", locus_status), grepl("MUT", locus_status))

    # Define a function to select and optionally normalize size_table_outer
    select_and_normalize <- function(fragment_obj, pattern, normalize = FALSE) {
      table <- as.data.frame(fragment_obj) %>% filter(grepl(pattern,
                                                     locus_status,
                                                     ignore.case = TRUE))
      if (normalize) {
        set.seed(123)
        table <- sample_n(table, nrow(size_table_inner))
      }
      table
    }

    # Apply conditions based on fragment_type
    if (fragment_type == "mut_vs_ref") {
      size_table_outer <- select_and_normalize(fragment_obj, "REF")
    } else if (fragment_type %in% c("mut_vs_outer",
                                    "mut_vs_ref_norm",
                                    "mut_vs_outer_norm")) {
      normalize <- grepl("norm", fragment_type) # normalisation check
      pattern <- sub(".*vs_([^_]*).*", "\\1", fragment_type) # check type
      size_table_outer <- select_and_normalize(fragment_obj,
                                               pattern,
                                               normalize)
    }

    summary_size_table_outer <- as.data.frame(size_table_outer %>%
                                                group_by(width) %>%
                                                summarize(count = n()))

    summary_size_table_inner <- as.data.frame(size_table_inner %>%
                                                group_by(width) %>%
                                                summarize(count = n()))

    summary_size_table_inner$MUTANT <- "true"
    summary_size_table_outer$MUTANT <- "false"

    size_merged_df <- merge(summary_size_table_outer,
                            summary_size_table_inner,
                            all = TRUE)

    colnames(size_merged_df) <- c("SIZE", "COUNT", "MUTANT")

    sizeCharacterisationSummary <- size_merged_df %>%
      mutate(SIZE.ROUNDED = plyr::round_any(SIZE, accuracy = 5L)) %>%
      group_by(MUTANT, SIZE.ROUNDED) %>%
      summarise(COUNT = sum(COUNT), .groups = "drop_last") %>%
      ungroup() %>%
      mutate(TOTAL = sum(COUNT)) %>%
      mutate(PROPORTION = COUNT / TOTAL)

    result <- sizeCharacterisationSummary %>%
      mutate(MUTANT_LABEL = as.factor(ifelse(MUTANT,
                                             "Mutation Fragment",
                                             "Reference Fragment")))

  }

  return(result)

}
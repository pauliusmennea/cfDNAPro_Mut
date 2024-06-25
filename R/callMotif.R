
#' Call Motif
#' @import dplyr tibble purrr tidyr
#' @param fragment_obj
#' @param genome_label
#' @param integrate_mut
#' @param motif_type
#' @param motif_length
#' @param ...
#'
#' @return tibble object
#' @export
#'
#' @examples
callMotif <- function(fragment_obj,
                      genome_label  = "hg19",
                      integrate_mut = FALSE,
                      motif_type = "s",
                      motif_length = 3L,
                      ...) {

  motif_length <- as.numeric(motif_length)

  message("Started to extract ", paste0(motif_type, motif_length), " motif...")
  message("You reference genome was set to ", genome_label, "!")

  # select the correct reference genome.

  if (genome_label == "hg19") {

    bsgenome_obj <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19

  }else if (genome_label == "hg38") {

    bsgenome_obj <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38

  }else if (genome_label == "hg38-NCBI") {

    bsgenome_obj <- BSgenome.Hsapiens.NCBI.GRCh38::BSgenome.Hsapiens.NCBI.GRCh38

  }

  if (integrate_mut == TRUE) {

    # Extract references and mutations based on locus_status
    gr_ref <- fragment_obj[grepl("REF", fragment_obj$locus_status)]
    gr_mut <- fragment_obj[grepl("MUT:concordant|MUT:single_read",
                                 fragment_obj$locus_status)]

    # Process motifs for reference and mutation fragments
    motif_calls_ref <- processMotif(fragment_obj = gr_ref,
                                    bsgenome_obj = bsgenome_obj,
                                    motif_type = motif_type,
                                    motif_length = motif_length)

    motif_calls_mut <- processMotif(fragment_obj = gr_mut,
                                    bsgenome_obj = bsgenome_obj,
                                    motif_type = motif_type,
                                    motif_length = motif_length)

    # Join reference and mutation motif calls
    result_frac <- full_join(motif_calls_ref, motif_calls_mut, by = "motif")

    # Rename with specific column names
    result_frac <- result_frac %>%
      dplyr::rename(
        motif = motif,
        n_ref = n.x,
        fraction_ref = fraction.x,
        n_mut = n.y,
        fraction_mut = fraction.y
      )

  } else {
    # Processing when integrate_mut is FALSE
    result_frac <- processMotif(fragment_obj, bsgenome_obj,
                                motif_type, motif_length)
  }

  message("Job completed successfully. ")
  return(result_frac)

}

# Helper function
processMotif <- function(fragment_obj, bsgenome_obj, motif_type, motif_length) {

  # Extract the motif to a vector
  motif <- get_motif(obj = fragment_obj,
                     genome = bsgenome_obj,
                     motif_type = motif_type,
                     motif_length = motif_length)

  # Summarize the vector and convert to tibble format
  result <- table(motif) %>%
    tibble::as_tibble()

  # Create a vector of elements
  pos <- c("C", "G", "A", "T")
  base_index <- seq.int(1, motif_length, by = 1)

  letter_list <- lapply(base_index, function(x, y) return(y), y = pos) %>%
    setNames(paste0("base", base_index))

  motif_ref <- purrr::cross_df(letter_list) %>%
    tidyr::unite(col = "motif", all_of(base_index), sep = "")

  # Report abnormal motifs.
  # Handle ambiguous motif(s)
  ambiguous_motif <- dplyr::anti_join(result, motif_ref, by = "motif")
  missing_motif <- dplyr::anti_join(motif_ref, result, by = "motif")

  if (nrow(ambiguous_motif) != 0) {
    message("ambiguous motif (i.e., 'N' base exists) detected: ")
    print(ambiguous_motif)

    result <- dplyr::filter(result, !stringr::str_detect(motif, "N"))
    message("Ambiguous motif(s) removed from result!")
  }

  # Handle missing motif(s)
  if (nrow(missing_motif) != 0) {
    message("Missing motif detected: ")
    print(missing_motif)

    result <- dplyr::right_join(result, motif_ref, by = "motif") %>%
      tidyr::replace_na(list(n = 0))

    message("Missing motif added back to the final result with count of 0!")
  }

  # Calculate the fraction
  result_frac <- result %>%
    dplyr::mutate(fraction = n / sum(n))

  # Set the factor level of motifs
  result_frac$motif <- factor(result_frac$motif, levels = sort(motif_ref$motif))

  return(result_frac)
}
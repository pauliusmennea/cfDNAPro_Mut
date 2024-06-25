#' Call trinucleotides and summarise the cfDNA information
#' for each target mutation locus.
#'
#' This function processes a GRanges object, summarizing cfDNA fragment
#' information for each target mutation locus.
#' It annotates each locus with the number and type of supporting fragments.
#' Each mismatch type is annotated with the median fragment length.
#' Consensus mismatch is determined by selecting the most frequent mismatch,
#' giving priority to the target mutation's ALT base over other bases.
#' Consensus mismatch is used to derive the trinucleotide substitution (SBS96).
#'
#' @import GenomicRanges
#' @import dplyr
#' @import BSgenome
#' @import IRanges
#' @import stringr
#'
#' @param gr GRanges object containing genomic ranges and associated data.
#'
#' @return dataframe with summarised mutational and trinucleotide data
#'
#' @examples
#' trinuc_df <-  callTrinucleotide(gr)
callTrinucleotide <- function(gr) {

  # Get ref genome info
  genome_info <- paste0("BSgenome.Hsapiens.UCSC.", unique(as.character(seqinfo(gr)@genome)))

  # Extract the BSgenome object from the package
  genome <- get(genome_info, envir = asNamespace(genome_info))

  # Convert GRanges object to dataframe and remove rows with NA in the target column
  gr_df <- as.data.frame(gr)
  gr_df <- gr_df[complete.cases(gr_df$target_mutation), ]

  # Create the new consensus_mismatch column
  gr_df <- gr_df %>%
    mutate(
      consensus_mismatch = paste0(
        sub("^((?:[^:]*:){3}).*", "\\1", locus_info),
        ifelse(
          grepl("REF", locus_status), "REF",
          ifelse(
            grepl("MUT:single_read|MUT:concordant", locus_status), "MUT",
            ifelse(
              grepl("other_base:concordant", locus_status),
              "other_base_concordant",
              ifelse(
                grepl("other_base:single_read", locus_status),
                "other_base_single_read",
                ifelse(
                  grepl("discordant", locus_status), "discordant", ""
                )
              )
            )
          )
        )
      )
    )

  # Function to process each row and conditionally remove the matching character
  process_discordant <- function(
      locus_status,
      consensus_mismatch,
      target_mutation) {
    if (grepl("discordant", locus_status)) {
      # Extract the character after the second colon in target_mutation
      char_in_target <- sub("^[^:]*:[^:]*:([^:]*):.*$", "\\1", target_mutation)
      # Extract the part of consensus_mismatch up to the second colon
      # and the part after the second colon
      prefix <- sub("^([^:]*:[^:]*:)[^:]*(:.*)$", "\\1", consensus_mismatch)
      suffix <- sub("^([^:]*:[^:]*:)[^:]*(:.*)$", "\\2", consensus_mismatch)
      # Extract the string between the second
      # and third colons in consensus_mismatch
      between_colons <- sub(
                            "^[^:]*:[^:]*:([^:]*):.*$",
                            "\\1",
                            consensus_mismatch)
      # Remove the character if it matches
      new_between_colons <- gsub(char_in_target, "", between_colons)
      # Reconstruct the consensus_mismatch value
      new_consensus_mismatch <- paste0(prefix, new_between_colons, suffix)
      return(new_consensus_mismatch)
    } else {
      return(consensus_mismatch)
    }
  }

  # Update the consensus_mismatch column conditionally row by row
  gr_df$consensus_mismatch <- mapply(
                                     process_discordant,
                                     gr_df$locus_status,
                                     gr_df$consensus_mismatch,
                                     gr_df$target_mutation,
                                     USE.NAMES = FALSE)

  # Function to randomly remove one character if
  # there are still two characters between colons
  process_bases <- function(
      consensus_mismatch,
      target_mutation) {
    # Extract the string between the second
    # and third colons in consensus_mismatch
    between_colons <- sub(
                          "^[^:]*:[^:]*:([^:]*):.*$",
                          "\\1",
                          consensus_mismatch)
    # Extract the character after the third colon in target_mutation
    char_after_third_colon <- sub(
                                  "^[^:]*:[^:]*:[^:]*:([^:]).*$",
                                  "\\1",
                                  target_mutation)
    # Check if there are two characters
    if (nchar(between_colons) == 2) {
      chars <- strsplit(between_colons, "")[[1]]
      # Check if either character matches the character
      # after the third colon in target_mutation
      if (char_after_third_colon %in% chars) {
        char_to_remove <- char_after_third_colon
      } else {
        # Randomly select one character to remove
        char_to_remove <- sample(chars, 1)
      }
      # Remove the selected character
      new_between_colons <- gsub(
                                 char_to_remove,
                                 "",
                                 between_colons, fixed = TRUE)
      # Reconstruct the consensus_mismatch value
      prefix <- sub("^([^:]*:[^:]*:)[^:]*(:.*)$", "\\1", consensus_mismatch)
      suffix <- sub("^([^:]*:[^:]*:)[^:]*(:.*)$", "\\2", consensus_mismatch)
      new_consensus_mismatch <- paste0(prefix, new_between_colons, suffix)
      return(new_consensus_mismatch)
    } else {
      return(consensus_mismatch)
    }
  }

  # Apply the removal function row by row to the discordant types
  # If one of the bases matches the ALT base of the target mutation
  # It will be prioritised over the other non-reference base
  gr_df$consensus_mismatch <- mapply(process_bases,
                          gr_df$consensus_mismatch,
                          gr_df$target_mutation,
                          USE.NAMES = FALSE)

  # Add variant fragment summary per mutation locus
  df1 <- gr_df %>%
    group_by(target_mutation) %>%
    summarize(
      CO_MUT = sum(locus_status %in% c("MUT:concordant")),
      SO_MUT = sum(locus_status %in% c("MUT:single_read")),
      CO_REF = sum(locus_status %in% c("REF:concordant")),
      SO_REF = sum(locus_status %in% c("REF:single_read")),
      DO = sum(locus_status %in% c("MUT:discordant")),
      SO_OTHER = sum(locus_status %in% c("other_base:single_read")),
      CO_OTHER = sum(locus_status %in% c("other_base:concordant"))
    )

  # Add fragment length summary per mutation locus
  df2 <- gr_df %>%
    group_by(target_mutation) %>%
    summarize(
      CO_MUT_flength = median(width[locus_status == "MUT:concordant"]),
      SO_MUT_flength = median(width[locus_status == "MUT:single_read"]),
      CO_REF_flength = median(width[locus_status == "REF:concordant"]),
      SO_REF_flength = median(width[locus_status == "REF:single_read"]),
      DO_flength = median(width[locus_status == "MUT:discordant"]),
      SO_OTHER_flength = median(width[locus_status == "other_base:single_read"]),
      CO_OTHER_flength = median(width[locus_status == "other_base:concordant"])
    )

  # Merge the two summaries
  merged_table_df <- merge(df1, df2, by = 'target_mutation')

  set.seed(123)  # For reproducibility in random selection

  # Mapping definitions
  mapping <- list(
    CO_MUT = "MUT",
    SO_MUT = "MUT",
    DO = "discordant",
    SO_OTHER = "other_base_single_read",
    CO_OTHER = "other_base_concordant"
  )

  # Function to get the highest value column(s)
  # and select one based on the priority rules
  get_highest_column <- function(row) {
    values <- row[c("CO_MUT", "SO_MUT", "DO", "SO_OTHER", "CO_OTHER")]
    max_value <- max(values)
    highest_columns <- names(values)[values == max_value]

    # If the choice is between SO_OTHER, CO_OTHER,
    # and DO with equal values, randomly choose one
    if (length(highest_columns) > 1) {
      selected_column <- sample(highest_columns, 1)
    } else {
      selected_column <- highest_columns
    }

    # Always prioritize CO_MUT or SO_MUT
    if ("CO_MUT" %in% highest_columns) {
      selected_columns <- "CO_MUT"
    }
    if ("SO_MUT" %in% highest_columns) {
      selected_columns <- "SO_MUT"
    }

    return(selected_column)
  }

  # Initialize the consensus_mismatch column in merged_table_df
  merged_table_df$consensus_mismatch <- NA

  # Iterate over each row in merged_table_df
  for (i in 1:nrow(merged_table_df)) {
    target_mutation <- merged_table_df$target_mutation[i]
    selected_column <- get_highest_column(merged_table_df[i, ])

    # Find matching rows in gr_df by target_mutation
    matching_rows <- gr_df[gr_df$target_mutation == target_mutation, ]

    # Define the required status based on the selected column
    required_status <- mapping[[selected_column]]

    # Filter the matching rows based on the required status
    if (required_status == "MUT") {
      filtered_rows <- matching_rows[
        grepl("MUT", matching_rows$locus_status), ]
    } else if (required_status == "discordant") {
      filtered_rows <- matching_rows[
        grepl("discordant", matching_rows$locus_status), ]
    } else if (required_status == "SO_OTHER") {
      filtered_rows <- matching_rows[
        grepl("other_base_single_read", matching_rows$locus_status),
      ]
    } else if (required_status == "CO_OTHER") {
      filtered_rows <- matching_rows[
        grepl("other_base_concordant", matching_rows$locus_status),
      ]
    }

    # Randomly select one of the filtered rows if multiple matches exist
    if (nrow(filtered_rows) > 1) {
      selected_row <- filtered_rows[sample(1:nrow(filtered_rows), 1), ]
    } else if (nrow(filtered_rows) == 1) {
      selected_row <- filtered_rows
    } else {
      selected_row <- NULL
    }

    # Update the consensus_mismatch column
    # in merged_table_df if a match is found
    if (!is.null(selected_row)) {
      merged_table_df$consensus_mismatch[i] <- selected_row$consensus_mismatch
    }
  }

  # Using complete.cases to focus on the consensus_mismatch column
  merged_table_df <- merged_table_df[complete.cases(merged_table_df$consensus_mismatch), ]


  # Split consensus_mismatch into three fields: chr, start, ALT Base
  split_values <- strsplit(merged_table_df$consensus_mismatch, ":")

  # Form a dataframe out of chr, start, end
  # (required for GRanges to get ref Bases)
  tri_df <- data.frame(
    chromosome = sapply(split_values, "[", 1),
    start = as.integer(sapply(split_values, "[", 2)) - 1,
    end = as.integer(sapply(split_values, "[", 2)) + 1
  )

  # Form a GRanges to obtain reference trinucleotide
  tri_loci_gr <- GRanges(seqnames = tri_df[, 1],
                         ranges = IRanges(
                                          start = tri_df[, 2],
                                          end = tri_df[, 3]))

  # Get Vector of trinucleotide
  ref_tri_bases <- as.vector(BSgenome::getSeq(genome, names = tri_loci_gr))

  # Add reference trinucleotide to table
  merged_table_df$ref_tri <- ref_tri_bases

  # Stratify into more columns to later obtain accurate trinucleotide mutations
  merged_table_df$mut_base <- sapply(split_values, "[", 3)
  merged_table_df$ref_base <- substr(merged_table_df$ref_tri, 2, 2)

  # Split df into two parts: A/G and C/T
  merged_table_df_1 <- merged_table_df %>% filter(ref_base %in% c("A", "G"))
  merged_table_df_2 <- merged_table_df %>% filter(ref_base %in% c("C", "T"))

  # Reverse complement for A/G
  merged_table_df_1$ref_tri <- chartr(
                                      "ATGC",
                                      "TACG",
                                      merged_table_df_1$ref_tri)
  merged_table_df_1$ref_base <- chartr(
                                       "ATGC",
                                       "TACG",
                                       merged_table_df_1$ref_base)
  merged_table_df_1$mut_base <- chartr(
                                       "ATGC",
                                       "TACG",
                                       merged_table_df_1$mut_base)

  merged_table_df_full <- rbind(merged_table_df_1, merged_table_df_2)

  merged_table_df_full$left <- substr(merged_table_df_full$ref_tri, 1, 1)
  merged_table_df_full$right <- substr(merged_table_df_full$ref_tri, 3, 3)

  merged_table_df_full$SBS96 <- paste0(
    merged_table_df_full$left,
    sep = "[",
    merged_table_df_full$ref_base,
    sep = ">",
    merged_table_df_full$mut_base,
    sep = "]",
    merged_table_df_full$right
  )

  # Add chr and pos columns
  merged_table_df_full$chr <- str_split(
                                        merged_table_df_full$consensus_mismatch,
                                        ":",
                                        simplify = TRUE)[, 1]
  merged_table_df_full$pos <- str_split(
                                        merged_table_df_full$consensus_mismatch,
                                        ":",
                                        simplify = TRUE)[, 2]

  # Remove unnecessary columns
  merged_table_df_full <- merged_table_df_full %>%
    select(-mut_base, -ref_base, -ref_tri, -left, -right, -chr, -pos)

  return(merged_table_df_full)
}

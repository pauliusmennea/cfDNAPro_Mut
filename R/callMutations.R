callMutations <- function(
                      bamfile,
                      use_names,
                      strand_mode,
                      galp_what,
                      galp_tag,
                      genome,
                      chromosome_to_keep,
                      galp_flag,
                      max_depth,
                      min_base_quality,
                      min_mapq,
                      min_nucleotide_depth,
                      min_minor_allele_depth,
                      pileupYieldSize,
                      frag) {
      # Define mutation coordinates for readGAlignements()
      bed_full <- pileupMismatches(
                            bamfile,
                            genome,
                            chromosome_to_keep,
                            galp_flag,
                            max_depth,
                            min_base_quality,
                            min_mapq,
                            min_nucleotide_depth,
                            min_minor_allele_depth,
                            pileupYieldSize)

      loci_df <- subset(bed_full, chr %in% chromosome_to_keep)

      target_muts <- bed_full %>% select(chr, start, ref, alt)

      target_muts$target <- paste(target_muts$chr,
                                  target_muts$start,
                                  target_muts$ref,
                                  target_muts$alt,
                                  sep = ":")

      target_muts$locus <- paste(target_muts$chr, target_muts$start, sep = ":")

      target_muts <- target_muts %>% select(locus, target)

      which_processed <- make_granges(loci = bed_full)

      mismatch_df <- processMismatches(
                      which_loci = which_processed,
                      loci_df = loci_df,
                      bamfile = bamfile,
                      galp_flag = galp_flag,
                      galp_mapqFilter = galp_mapqFilter,
                      chromosome_to_keep = chromosome_to_keep)

    # Match the values and if they match by coordinate
    # Add target_muts$target as target column 'target_mutation' in mismatch_df
    mismatch_df$locus <- sub("^([^:]+:[^:]+):.*", "\\1", mismatch_df$locus_info)

    target_muts <- loci_df %>% select(chr, start, ref, alt)

    target_muts$target <- paste(
                                target_muts$chr,
                                target_muts$start,
                                target_muts$ref,
                                target_muts$alt,
                                sep = ":")

    target_muts$locus <- paste(target_muts$chr, target_muts$start, sep = ":")

    target_muts <- target_muts %>% select(locus, target)

    # Match values and add 'target_mutation' column to the first dataframe
    mismatch_df <- merge(mismatch_df, target_muts, by = "locus", all.x = TRUE)

  # Tidy up the reference base when discordant
  # Function to replace the specific R in locus_info based on ref_base
  replace_specific_R <- function(locus_status, locus_info, target) {
    if (locus_status == "MUT:discordant") {
      ref_base <- sub("^[^:]*:[^:]*:([^:]).*", "\\1", target)
      sub("R", ref_base, locus_info)
    } else {
      locus_info
    }
  }

  # Update the locus_info column conditionally row by row
  mismatch_df$locus_info <- mapply(replace_specific_R,
                                 mismatch_df$locus_status,
                                 mismatch_df$locus_info,
                                 mismatch_df$target,
                                 USE.NAMES = FALSE)

    # Remove .1, .2, etc., annotations from fragment_id column
    mismatch_df$fragment_id <- gsub("\\.\\d+$", "", mismatch_df$fragment_id)

    # GRanges to dataframe
    # This adds .1 to row names whenever there is a duplicate fragment ID
    metadata <- as.data.frame(names(frag))

    colnames(metadata) <- "fragment_event"

    metadata$fragment_event <- make.unique(metadata$fragment_event, sep = ".")

    # Allow duplicate row names by adding .1 .2 .3
    fragment_gr_df <- as.data.frame(frag, use.outer.mcols = TRUE,
                                    row.names = metadata$fragment_event)

    fragment_gr_df$fragment_id <- rownames(fragment_gr_df)

    # Remove suffixes like .1, .2, .3, etc.
    fragment_gr_df$fragment_id <- sub("\\.\\d+$", "", fragment_gr_df$fragment_id)

    # Remove duplicates in the first dataframe - new solution
    fragment_gr_df <- fragment_gr_df %>% distinct(fragment_id, .keep_all = TRUE)

    # Merge dataframes using base R merge - new solution
    merged_df <- suppressWarnings(merge(
      fragment_gr_df,
      mismatch_df,
      by = "fragment_id",
      all.x = TRUE,
      all.y = FALSE
    ))

    # Tidy up dataframe:
    # Extract the first letter following the first ':' in the target column
    replacement_chars <- sapply(strsplit(merged_df$target, ":"),
                                function(x) substr(x[3], 1, 1))

    # Replace 'REF' in the locus_info column with the extracted character
    merged_df$locus_info <- mapply(function(info, replacement) {
      sub("REF", replacement, info)
    }, merged_df$locus_info, replacement_chars)

    # Rename target column to target_mutation
    names(merged_df)[names(merged_df) == "target"] <- "target_mutation"

    # Convert df to GRanges:
    # Extract Seqinfo from the existing frag GRanges object and add it back
    seqinfo_frag <- seqinfo(frag)
    merged_gr <- suppressWarnings(makeGRangesFromDataFrame(
                                                merged_df,
                                                keep.extra.columns = TRUE,
                                                seqinfo = seqinfo_frag))

    gr_meta <- suppressWarnings(mcols(merged_gr))

    gr_meta$locus_info[is.na(gr_meta$locus_info)] <- "outer_fragment"
    gr_meta$locus_status[is.na(gr_meta$locus_status)] <- "outer_fragment"
    mcols(merged_gr) <- gr_meta

    return(merged_gr)

  }
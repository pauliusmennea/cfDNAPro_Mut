#' Internal Functions for processing mutational fragment-level information
#'
#'

check_mutfile_columns <- function(df) {

  # Check if the data frame has four columns
  if (ncol(df) != 4) {
    return("Data frame does not have four columns")
  }

  # Check if the column names are correct
  col_names <- colnames(df)
  if (!all(col_names %in% c("chr", "pos", "ref", "alt"))) {
    return("Column names are incorrect")
  }

  # Check if the column types are correct
  col_types <- sapply(df, class)
  if (!all(col_types == c("character", "integer", "character", "character"))) {
    return("Column types are incorrect")
  }

  # If all checks pass, return "OK"
  return("Mutation file OK!")
}



#' Read in the .tsv mutation file
#'
#' @importFrom utils read.table


read_mutation_file <- function(mutation_file) {

  # Check if file header is present
  first_line <- readLines(mutation_file, n = 1)
  has_header <- first_line == paste0(c("chr", "pos", "ref", "alt"),
                                     collapse = "\t")

  # Read in the tab separated mutation file
  message("Reading in the provided mutation file ...")

  if (has_header) {
    bed <- read.table(mutation_file,
                      header = TRUE,
                      sep = "\t",
                      stringsAsFactors = FALSE,
                      quote = "")
  } else {
    bed <- read.table(mutation_file,
                      header = FALSE,
                      sep = "\t",
                      stringsAsFactors = FALSE, quote = "",
                      col.names = c("chr", "pos", "ref", "alt"))
  }

  # Check whether the file format is correct
  result <- check_mutfile_columns(bed)

  if (result == "Mutation file OK!") {
    message("Mutation file OK!")
    # Process the file
    bed["start"] <- bed["pos"]
    bed["end"] <- bed["pos"]
    bed_full <- bed %>% select(chr, start, end, ref, alt)

    return(bed_full)

    # Error message
  } else {
    message("Incorrect Mutation File Format")
    message("Please provide a tab separated file with 4 columns:")
    message("chr	pos	ref	alt")
    message("Column header is not required")
  }
}


#' Read and Process the provided Mutation File
#'
#' @importFrom gsubfn strapply
#' @importFrom purrr simplify
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom magrittr %>%
#'
#' @param bamfile
#' @param use_names
#' @param chromosome_to_keep
#' @param strand_mode
#' @param genome_label
#' @param galp_flag
#' @param galp_what
#' @param galp_tag
#' @param galp_mapqFilter
#' @param mutation_file
#' @param frag
#'
#' @return Returns a GRanges object with curated mutational information
#' @export
#'
#' @examples
readMutations <- function(
                      bamfile,
                      use_names,
                      chromosome_to_keep,
                      strand_mode,
                      genome_label,
                      galp_flag,
                      galp_what,
                      galp_tag,
                      galp_mapqFilter,
                      mutation_file,
                      frag) {
  # Define mutation coordinates for readGAlignements()
  loci_df <- read_mutation_file(mutation_file)

  loci_df <- subset(loci_df, chr %in% chromosome_to_keep)

  which_loci <- make_granges(loci = loci_df)

  mismatch_df <- processMismatches(
                      which_loci = which_loci,
                      loci_df = loci_df,
                      bamfile = bamfile,
                      galp_flag = galp_flag,
                      galp_mapqFilter = galp_mapqFilter,
                      chromosome_to_keep = chromosome_to_keep)


  # Match the values, and if they match by coordinate
  # add target_muts$target as target column 'target_mutation' in mismatch_df
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

  # allow duplicate row names by adding .1 .2 .3
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


#' Subsetting alignments to mutational positions
#' @import magrittr
#' @import GenomeInfoDb
#' @import GenomicAlignments
#' @import S4Vectors
#' @import Rsamtools

bam_to_galp_mut <- function(bamfile,
                         use_names = TRUE,
                         param = galp_param,
                         chromosome_to_keep = FALSE,
                         strand_mode = 1,
                         genome = NA_character_) {
  # Check parameters
  stopifnot(file.exists(bamfile))
  stopifnot(isSingleStringOrNA(genome) || is(genome, "Seqinfo"))

  # Read bam into galp
  message("Reading bam into galp...")

  galp <- readGAlignmentPairs(file = bamfile,
                              use.names = use_names,
                              strandMode = strand_mode,
                              param = param,
                              with.which_label = TRUE)
  # add genome information
  if (isSingleStringOrNA(genome)) {
    genome <- Seqinfo(genome = genome)
  }
  seqinfo(galp) <- merge(seqinfo(galp), genome)

  # strandMode should be one for downstream operations
  stopifnot(GenomicAlignments::strandMode(galp) == 1)

  # only keep needed seqnames
  if (!isFALSE(chromosome_to_keep)) {
    galp <- keepSeqlevels(galp, chromosome_to_keep, pruning.mode = "coarse")

  }

  message("Curating seqnames and strand information...")
  # remove read pairs without correct seqnames and strand information
  galp2 <- galp[!is.na(GenomicAlignments::seqnames(galp))]
  galp3 <- galp2[GenomicAlignments::strand(galp2) != "*"]

  return(galp3)

}

#' Generate GRanges for specific positions
#' @importFrom IRanges IRanges

make_granges <- function(loci) {
  loci_gr <- GenomicRanges::GRanges(seqnames = loci[, 1],
                     ranges = IRanges::IRanges(start = loci[, 2],
                                      end = loci[, 3]))
  return(loci_gr)
}

#' Process and curate mutational fragment-level information
process_mutation_fragments <- function(
    bamfile = bamfile,
    use_names = use_names,
    chromosome_to_keep = chromosome_to_keep,
    strand_mode = strand_mode,
    genome = genome,
    galp_flag = galp_flag,
    galp_what = galp_what,
    galp_tag = galp_tag,
    galp_mapqFilter = galp_mapqFilter,
    mutation_file = mutation_file,
    call_mutations = call_mutations,
    pileup_max_depth = pileup_max_depth,
    pileup_min_base_quality = pileup_min_base_quality,
    pileup_min_mapq = pileup_min_mapq,
    pileup_min_nucleotide_depth = pileup_min_nucleotide_depth,
    pileup_min_minor_allele_depth = pileup_min_minor_allele_depth,
    pileupYieldSize = pileupYieldSize,
    frag = frag) {

  # Check if mutation file is provided and use it for annotation
  if (!is.null(mutation_file)) {
    frag <- suppressWarnings(readMutations(
                          bamfile = bamfile,
                          use_names = use_names,
                          chromosome_to_keep = chromosome_to_keep,
                          strand_mode = strand_mode,
                          genome = genome,
                          galp_flag = galp_flag,
                          galp_what = galp_what,
                          galp_tag = galp_tag,
                          galp_mapqFilter = galp_mapqFilter,
                          mutation_file = mutation_file,
                          frag = frag))

  } else if (call_mutations == TRUE && is.null(mutation_file)) {

    # Call mutations if no mutation file is provided and call_mutations is TRUE
    frag <- callMutations(
                    bamfile = bamfile,
                    use_names = use_names,
                    strand_mode = strand_mode,
                    galp_what = galp_what,
                    galp_tag = galp_tag,
                    genome = genome,
                    chromosome_to_keep = chromosome_to_keep,
                    galp_flag = galp_flag,
                    max_depth = pileup_max_depth,
                    min_base_quality = pileup_min_base_quality,
                    min_mapq = galp_mapqFilter,
                    min_nucleotide_depth = pileup_min_nucleotide_depth,
                    min_minor_allele_depth = pileup_min_minor_allele_depth,
                    yieldSize = pileupYieldSize,
                    frag = frag)

  } else if (call_mutations == FALSE && is.null(mutation_file)) {
    frag <- frag
  }
  return(frag)
}

#' Function to process BAM files based on mutational or general alignment data
#' @import magrittr
#' @import GenomeInfoDb
#' @import GenomicAlignments
#' @import S4Vectors
#' @import Rsamtools

process_bam_file <- function(
                             bamfile = bamfile,
                             mutation_file = mutation_file,
                             mut_fragments_only = mut_fragments_only,
                             use_names = use_names,
                             chromosome_to_keep = chromosome_to_keep,
                             strand_mode = strand_mode,
                             genome_name = genome_name,
                             call_mutations = call_mutations,
                             galp_flag = galp_flag,
                             galp_what = galp_what,
                             galp_tag = galp_tag,
                             galp_mapqFilter = galp_mapqFilter) {

  # Check if only mutations are to be targeted
  if (mut_fragments_only) {

    # Process mutation-specific alignments
    loci_df <- read_mutation_file(mutation_file)
    loci_df <- subset(loci_df, chr %in% chromosome_to_keep)
    which_loci <- make_granges(loci = loci_df)

    # Initialize general parameters
    galp_param <- Rsamtools::ScanBamParam(
      flag = galp_flag,
      what = galp_what,
      tag = galp_tag,
      mapqFilter = galp_mapqFilter,
      which = which_loci)

    # Fetch mutation specific alignments
    galp <- bam_to_galp_mut(
      bamfile = bamfile,
      use_names = use_names,
      chromosome_to_keep = chromosome_to_keep,
      strand_mode = strand_mode,
      genome = genome_name,
      param = galp_param)

  } else {

    # Initialize general parameters
    galp_param <- Rsamtools::ScanBamParam(
      flag = galp_flag,
      what = galp_what,
      tag = galp_tag,
      mapqFilter = galp_mapqFilter)

    # Process all alignments
    galp <- bam_to_galp2(
      bamfile = bamfile,
      use_names = use_names,
      chromosome_to_keep = chromosome_to_keep,
      strand_mode = strand_mode,
      genome = genome_name,
      param = galp_param)
  }

  # Optionally remove alignments with INDELs
  if (!is.null(mutation_file) || call_mutations == "TRUE") {

    # Identify and exclude indel fragments
    indel_aln <- grepl("[ID]", cigar(galp@first)) |
      grepl("[ID]", cigar(galp@last))
    galp <- galp[!indel_aln]
  }

  # Return the processed data
  return(galp)
}

#' Remove soft-clipped bases from read sequences
#' @importFrom gsubfn strapply
#' @importFrom purrr simplify
#'
#' @param vdf
#'
#' @return a dataframe
#' @export
#'
#' @examples
clipReadSeq <- function(vdf = NULL) {

    # set up extra columns with null values
    vdf$softclip <- rep(0, nrow(vdf))
    vdf$leftclip <- rep(0, nrow(vdf))
    vdf$rightclip <- rep(0, nrow(vdf))
    vdf$clipseq   <- as.vector(vdf$seq)

    # we take a copy of the cigar string
    # and remove hard clipping information
    # for convenience
    vdf$cigar2 <- gsub("\\d+H", "", vdf$cigar)

    # Record the column indices for cigar columns
    cidx <- grep("cigar", colnames(vdf))[1]
    cidx2 <- grep("cigar", colnames(vdf))[2]

    # get the subset of rows with soft-clipping
    clipr <- grep("S", vdf$cigar)
    svdf <- vdf[clipr, ]

    # calculate the number of bases soft clipped from
    # the left of the read by extracting the number
    # adjacent to any S at the beginning of the string
    svdf$leftclip <- apply(svdf, 1, function(x) {
      sum(as.numeric(strapply(as.character(x[cidx2]), "^(\\d+)S", simplify = c)))
    })

    # calculate the number of bases soft clipped from
    # the right of the read by extracting the number
    # adjacent to any S at the end of the string
    svdf$rightclip <- apply(svdf, 1, function(x) {
      sum(as.numeric(strapply(as.character(x[cidx2]), "(\\d+)S$", simplify = c)))
    })

    # extract the sequence in the middle of leftclip and rightclip
    # this is the sequence that has actually aligned
    svdf$clipseq <- substr(svdf$seq,
                        as.numeric(svdf$leftclip) + 1,
                        nchar(as.vector(svdf$seq)) - as.numeric(svdf$rightclip))

    # substitute in the values
    # vdf$softclip[clipr] <- svdf$softclip
    vdf$leftclip[clipr] <- svdf$leftclip
    vdf$rightclip[clipr] <- svdf$rightclip
    vdf$clipseq[clipr] <- svdf$clipseq

    # calculate the length of the clipped sequence
    vdf$cliplen <- nchar(vdf$clipseq)

    # remove the second cigar string
    # we calculated above
    vdf <- vdf[, -cidx2]

    # return the data
    return(vdf)

}

#' Obtain internal position of mismatch within read
#'
#' @param read_id
#' @param chr_start
#' @param start_pos
#' @param read_strand
#' @param nm_tag
#' @param md_tag
#' @param cigar_string
#' @param clip_seq
#' @param which_locus
#' 
#' @return a dataframe
#' @export
#'
#' @examples

mismatchPosition <- function(
                    read_id,
                    chr_start,
                    start_pos,
                    read_strand,
                    nm_tag,
                    md_tag,
                    cigar_string,
                    clip_seq,
                    which_locus) {

    # First crucial step is to split the md string at
    # each substitution/deletion operation
    md_gsub <- gsub("([\\^]*[ACGT]+)[0]*", " \\1 ", md_tag)

    # Split each operation using strsplit
    md_spl  <- strsplit(md_gsub, "[ ]+")[[1]]
    this    <- as.integer()

    # since its a relatively time-consuming operation
    # use NM flag to calculate position of mismatches
    # only if there are mismatches
    if (nm_tag != 0) {
    this <- lapply(md_spl, function(y) {
        if (!is.na(as.numeric(y))) {
            o <- rep("M", as.numeric(y))
        } else if(length(grep("\\^", y)) > 0) {
            o <- rep("D", nchar(y) - 1)
        } else if (nchar(y) == 1) {
            o <- rep("MM", 1)
        }
    })

        # we get a list of consecutive M, D or MM occurences
        # separated each time a different type appeared
        # we combine the list of those vectors into one
        this <- do.call(c, this)

        # after this step, we have a vector of length =
        # read length and each position telling if it's a
        # match(M), mismatch(MM) or deletion (D)

        # next we find where in the read, the MM occurs
        this <- which(this == "MM")
    }

    # if there are no mismatches, we assign 0 to the read
    # 0 allows use to later differentiate REF reads
    if (length(this) == 0) {

    this <- 0

    }

    mismatch_base <- lapply(this, function(x) substr(clip_seq, x, x))

    genomic_coordinate <- paste0(
                                read_id,
                                sep = ":",
                                read_strand,
                                sep = "_",
                                chr_start,
                                sep = ":",
                                this + start_pos - 1,
                                sep = ":",
                                mismatch_base,
                                sep = ":read.position=",
                                this,
                                sep = ":",
                                which_locus)
  genomic_coordinate
}

#' Process Mismatch data
#' @importFrom gsubfn strapply
#' @import dplyr
#' @importFrom purrr simplify
#' @import Parallel
#' @import magrittr
#' @import GenomeInfoDb
#' @import GenomicAlignments
#' @import S4Vectors
#' @import Rsamtools
#'
#' @param which_loci
#' @param loci_df
#' @param bamfile
#' @param galp_flag
#' @param galp_mapqFilter
#'
#' @return a Dataframe of processes mismatches
#' @export
#'
#' @examples

processMismatches <- function(which_loci,
                              loci_df,
                              bamfile,
                              galp_flag,
                              galp_mapqFilter,
                              chromosome_to_keep) {

  ############################################################################
  # Process bam file paired-end reads and obtain mismatch information
  ############################################################################
  what <- c("seq", "mapq")

  param <- ScanBamParam(
                        which = which_loci,
                        flag = galp_flag,
                        what = what,
                        mapqFilter = galp_mapqFilter,
                        tag = c("NM", "MD"))

  galp_modified_raw <- readGAlignments(file = bamfile,
                                       use.names = TRUE,
                                       with.which_label = TRUE,
                                       param = param)

  # only keep needed seqnames
  if (!isFALSE(chromosome_to_keep)) {
    galp_modified_raw <- keepSeqlevels(galp_modified_raw,
                                       chromosome_to_keep,
                                       pruning.mode = "coarse")}
  message("Mismatch processing stage initiated ...")

  ############################################################################
  # Convert GRanges to DF and remove indels
  ############################################################################
  gr_df <- as.data.frame(galp_modified_raw, use.outer.mcols = TRUE)

  gr_df <- gr_df[!grepl("I", gr_df$cigar), ]

  gr_df <- gr_df[!grepl("D", gr_df$cigar), ]

  ############################################################################
  # Clip the read sequences to remove soft clipped bases to match the MD tags
  ############################################################################
  gr_df_clipped <- clipReadSeq(gr_df)
  message("Read sequences have been adjusted for soft-clipping... ")

  # Reads with .2, .3 are not necessarily overlapping multiple mutations
  # Instead they may simply by paired-end reads overlapping one mutation
  # We fix this by applying the read_id processing steps below
  gr_df_clipped$read_id <- rownames(gr_df_clipped)
  gr_df_clipped$read_id  <- gsub("\\.\\d+", "", gr_df_clipped$read_id)

  # Asigning strand information to each read id
  gr_df_clipped <- gr_df_clipped %>%
    mutate(read_id = paste(read_id, strand, sep = ""))

  # Adding a suffix to reads overlapping more than one locus
  gr_df_clipped$read_id <- make.unique(gr_df_clipped$read_id, sep = ".")

  # Removing unnecessary strand information
  # Now the suffix of the reads represents their overlap of multiple mutations
  # For exampe, If the read overlaps 2 mutations, then it will have .1 suffix.
  gr_df_clipped$read_id <- gsub("[+-]", "", gr_df_clipped$read_id)

  ############################################################################
  # Tidy up the clipped read dataframe and remove duplicate .1 annotation
  ############################################################################
  gr_df_clipped_reads <-  gr_df_clipped %>% select(
                          read_id,
                          seqnames,
                          start,
                          strand,
                          NM,
                          MD,
                          cigar,
                          clipseq,
                          which_label)

  # When fragments overlap multiple mutations, there will be >= 2 duplicates
  # Previously this was just .1 for gsub
  gr_df_clipped_reads$read_id <- gsub("\\.1",
                                      ".1",
                                      as.character(gr_df_clipped_reads$read_id))

  #gr_df_clipped_reads$read_id <- gsub("\\.[2-3]",
  #                                    ".2",
  #                                    as.character(gr_df_clipped_reads$read_id))
  #
  #gr_df_clipped_reads$read_id <- gsub("\\.[4-5]",
  #                                    ".3",
  #                                    as.character(gr_df_clipped_reads$read_id))
  #
  #gr_df_clipped_reads$read_id <- gsub("\\.[6-7]",
  #                                    ".4",
  #                                    as.character(gr_df_clipped_reads$read_id))
  #
  colnames(gr_df_clipped_reads) <- c("read_id",
                                     "chr_start",
                                     "start_pos",
                                     "read_strand",
                                     "nm_tag",
                                     "md_tag",
                                     "cigar_string",
                                     "clip_seq",
                                     "which_locus")

  ############################################################################
  # Annotated each read with mismatch information and generate a dataframe
  ############################################################################
  suppressWarnings(
    output_f <- purrr::pmap(gr_df_clipped_reads, mismatchPosition))

  message("Processing mismatch positions...")

  output_f_df <- as.data.frame(stri_list2matrix(output_f, byrow = TRUE))

  ############################################################################
  # Remove read names. They will be added back later
  ############################################################################
  output_f_df_clean <- as.data.frame(lapply(
                        output_f_df, function(x) gsub(".*\\_", "\\1", x)))

  ############################################################################
  # Remove unnecessary locus positions from mismatches
  ############################################################################
  output_f_df_clean <- as.data.frame(lapply(
                        output_f_df_clean, function(x) gsub("\\:c.*", "\\1", x)))

  ############################################################################
  # Add back read's locus of interest as a separate column
  ############################################################################
  output_f_df_clean$which_locus <- gsub(".*\\:c", "\\1c", output_f_df$V1)

  ############################################################################
  # Add back read ID with strand information (+/-) to a separate column
  ############################################################################
  output_f_df_clean$read_id <- gsub("\\_.*", "\\1", output_f_df$V1)

  # Get just the locus
  output_f_df_clean$locus <- gsub("-\\d+", "", output_f_df_clean$which_locus)

  # Define the V column names
  v_columns <- grep("^V\\d+$", names(output_f_df_clean), value = TRUE)

  # Define the number of rows per chunk
  chunk_size <- 10000

  # Calculate the number of chunks needed
  num_chunks <- ceiling(nrow(output_f_df_clean) / chunk_size)

  # Split the dataframe into chunks
  df_chunks <- split(
    output_f_df_clean, rep(
      1:num_chunks, each = chunk_size, length.out = nrow(output_f_df_clean)))

  # Function to process each chunk in parallel
  process_chunk <- function(chunk) {

    # Iterate over each row of the chunk
    for (i in 1:nrow(chunk)) {

      # Get the 'which_locus' value for the current row
      target_locus <- chunk$locus[i]

      # Check if V1 already contains the string "read.position=0"
      if (!grepl("read.position=0", chunk$V1[i], fixed = TRUE)) {

        # Check which V columns partially match with the 'which_locus' value
        matching_columns <- v_columns[
          sapply(v_columns, function(column) {
            grepl(target_locus, chunk[i, column], fixed = TRUE)
          })
        ]

        # If matching columns are found
        if (length(matching_columns) > 0) {

          # Move the first matching V column to V1
          chunk$V1[i] <- chunk[i, matching_columns[1]]

          # Set all other V columns (excluding V1) to NA
          chunk[i, v_columns[-1]] <- NA
        } else {

          # If no matching column is found,
          # set V1 as the 'which_locus' value plus "::read.position=0"
          chunk$V1[i] <- paste0(target_locus, "::read.position=0")

          # Set all V columns (excluding V1) to NA
          chunk[i, v_columns[-1]] <- NA
        }
      }
    }

    # Return the processed chunk
    return(chunk)
  }

  # Apply the function to each chunk in parallel
  cl <- makeCluster(4)  # Create a cluster with 4 cores

  # Export the required objects to cluster
  clusterExport(cl, c("v_columns"), envir = environment())

  processed_chunks <- parLapply(cl, df_chunks, process_chunk)

  stopCluster(cl)  # Stop the cluster

  # Merge the processed chunks back into a single dataframe
  output_f_df_clean <- do.call(rbind, processed_chunks)

  ############################################################################
  # Convert wide format into long format. This gets rid of the NAs.
  ############################################################################
  output_f_df_clean_gathered <- output_f_df_clean %>%
                                    gather(key = "key",
                                           value = "mismatch",
                                           starts_with("V"),
                                           na.rm = TRUE)

  ############################################################################
  # Tidy up the dataframe by selecting columns and removing unnecessary info
  ############################################################################
  output_f_df_clean_gathered <- output_f_df_clean_gathered %>%
                                    select(which_locus,
                                           read_id,
                                           mismatch)

  output_f_df_clean_gathered$which_locus <- gsub(
                                    "\\-.*",
                                    "\\1",
                                    output_f_df_clean_gathered$which_locus)

  ############################################################################
  # read.position=0 indicates that there there were no mismatches in a read
  ############################################################################

  ############################################################################
  # Thus, reads with read.position=0 read get assigned REF as their base
  ############################################################################

  ############################################################################
  # This suggests the read overlapped the locus of interest but had REF as base
  ############################################################################
  output_f_df_clean_gathered$allele <- ifelse(
                                        grepl(
                                        "read.position=0",
                                        output_f_df_clean_gathered$mismatch),
                                        "REF",
                                        "ALT")

  ############################################################################
  # Tidy up the dataframe by splitting the mismatch column into extra columns
  ############################################################################

  ############################################################################
  # Mismatch position within the read placed in its own column
  # Mismatch positions are in the 5'-> 3' direction regardless of read strand
  ############################################################################
  output_f_df_clean_gathered$read_position <- gsub(
                                        ".*\\:r",
                                        "\\1r",
                                        output_f_df_clean_gathered$mismatch)

  ############################################################################
  # Removing mismatch position information from the mismatch column
  ############################################################################
  output_f_df_clean_gathered$mismatch <- gsub(
                                        "\\:read.position.*",
                                        "\\1",
                                        output_f_df_clean_gathered$mismatch)

  ############################################################################
  # Listing SNV as REF if there were no mismatches found in the read
  ############################################################################
  output_f_df_clean_gathered$snv <- ifelse(grepl(
            "read.position=0",
            output_f_df_clean_gathered$read_position),
            paste0(output_f_df_clean_gathered$which_locus, sep = ":", "REF"),
            output_f_df_clean_gathered$mismatch)

  ############################################################################
  # Using the Mutation file to query the mismatch dataframe
  ############################################################################
  loci_df1 <- loci_df %>% select(chr, start, alt)

  ############################################################################
  # Obtain a vector of targeted mutations
  ############################################################################
  loci_vector <- apply(loci_df1, 1, function(x) paste(x, collapse = ":"))

  loci_vector <- gsub(" ", "", loci_vector)

  ############################################################################
  # Query whether each mismatch in every read matches the targeted mutation
  ############################################################################
  output_f_df_clean_gathered$match_mutation_bed <- ifelse(
                        output_f_df_clean_gathered$snv %in% loci_vector,
                        TRUE,
                        FALSE)

  ############################################################################
  # Assign fragment ID to each mismatch
  ############################################################################
  output_f_df_clean_gathered$fragment_id <- gsub(
                                    ":([^:]+)$",
                                    "",
                                    output_f_df_clean_gathered$read_id)

  ############################################################################
  # Place mismatch genomic positions into a separate column
  ############################################################################
  output_f_df_clean_gathered$snv_pos <- gsub(
                ":([^:]+)$",
                "",
                output_f_df_clean_gathered$snv)

  ############################################################################
  # Place read strand information into a separate column
  ############################################################################
  output_f_df_clean_gathered$strand <- gsub(
                    ".*\\:",
                    "",
                    output_f_df_clean_gathered$read_id)

  ############################################################################
  # Add SNV+Strand column to indicate which paired-end read had the mismatch
  ############################################################################
  output_f_df_clean_gathered$snv_read_specific <- paste0(
                    output_f_df_clean_gathered$snv,
                    sep = ":", output_f_df_clean_gathered$strand)

  ############################################################################
  # Obtain mutation of interest positions as a vector
  ############################################################################
  loci_df2 <- loci_df %>% select(chr, start)

  bed_vector_pos <- apply(
                          loci_df2,
                          1,
                          function(x) paste(x, collapse = ":"))

  bed_vector_pos <- gsub(" ", "", bed_vector_pos)

  ############################################################################
  # Focus on the mutation file positions to extract relevant information
  ############################################################################
  `%nin%` <- Negate(`%in%`)

  ############################################################################
  # Reads with mismatch positions absent in mutation file, get assigned REF
  ############################################################################
  output_f_df_clean_gathered$snv_read_specific <- ifelse(
    output_f_df_clean_gathered$snv_pos %nin% bed_vector_pos,
    paste0(output_f_df_clean_gathered$which_locus,
           sep = ":REF:",
           output_f_df_clean_gathered$strand),
    paste0(output_f_df_clean_gathered$snv,
           sep = ":",
           output_f_df_clean_gathered$strand))

  ############################################################################
  # New column to keep reads that match the positions of mutation file
  ############################################################################
  output_f_df_clean_gathered$match_mutation_bed_pos <- ifelse(
        output_f_df_clean_gathered$snv_pos %in% bed_vector_pos,
        TRUE,
        FALSE)

  ############################################################################
  # Dataframe tidy up
  ############################################################################
  final_output <- output_f_df_clean_gathered %>% select(
                      fragment_id,
                      snv_read_specific,
                      match_mutation_bed,
                      which_locus)

  final_output$snv_read_specific <- paste(final_output$snv_read_specific,
                                          sep = "|",
                                          final_output$match_mutation_bed)

  ############################################################################
  # Wide format to convert DF from read-level to fragment-level
  ############################################################################

  # Add a row number column to handle duplicates
  df <- final_output[, -3]

  # Convert fragment_id to character type
  df <- df %>%
    mutate(fragment_id = as.character(fragment_id))

  df$new_id <- paste(
    sub("\\..*$", "", df$fragment_id), ":", df$which_locus, sep = "")

  df_wide <- df %>%
          group_by(new_id) %>%
          mutate(row_num = row_number()) %>%
          ungroup() %>%
          pivot_wider(
            id_cols = new_id,
            names_from = row_num,
            values_from = snv_read_specific,
            names_prefix = "snv_read_specific_"
          ) %>%
          unnest(c(starts_with("snv_read_specific_")))

  df_wide$fragment_id <- gsub(":chr[0-9]+:[0-9]+", "", df_wide$new_id)

  df_wide$fragment_id <- make.unique(df_wide$fragment_id, sep = ".")

  final_output_wide <- df_wide %>% select(
    fragment_id,
    snv_read_specific_1,
    snv_read_specific_2)

  ############################################################################
  # Add fragment level information the mismatch and paired end read support
  ############################################################################
  final_output_wide$concordance <- ifelse(
    is.na(final_output_wide$snv_read_specific_1) |
    is.na(final_output_wide$snv_read_specific_2),
    NA,
    str_extract(
        final_output_wide$snv_read_specific_1, ".*:") == str_extract(
            final_output_wide$snv_read_specific_2, ".*:")
  )

  final_output_wide$bed_locus <- ifelse(
    !is.na(final_output_wide$snv_read_specific_1),
    gsub("(:[^:]+){2}$", "", final_output_wide$snv_read_specific_1),
    ifelse(
        !is.na(final_output_wide$snv_read_specific_2),
        gsub("(:[^:]+){2}$", "", final_output_wide$snv_read_specific_2),
    NA)
  )

  ############################################################################
  # Discern between: base = MUT file, base = REF, discordant R1/R2, other base
  ############################################################################
  final_output_wide$match  <- ifelse(
        grepl("TRUE", final_output_wide$snv_read_specific_1) &
        is.na(final_output_wide$snv_read_specific_2),
        gsub(":([^:]+)$", "", final_output_wide$snv_read_specific_1),
    ifelse(
        grepl("TRUE", final_output_wide$snv_read_specific_2) &
        is.na(final_output_wide$snv_read_specific_1),
        gsub(":([^:]+)$", "", final_output_wide$snv_read_specific_2),
    ifelse(
        grepl("REF", final_output_wide$snv_read_specific_1) &
        is.na(final_output_wide$snv_read_specific_2),
        paste0("REF"),
    ifelse(
        grepl("REF", final_output_wide$snv_read_specific_2) &
        is.na(final_output_wide$snv_read_specific_1),
        paste0("REF"),
    ifelse(
        grepl("TRUE", final_output_wide$snv_read_specific_1) &
        grepl("REF", final_output_wide$snv_read_specific_2),
        paste0("discordant"),
    ifelse(
        grepl("TRUE", final_output_wide$snv_read_specific_2) &
        grepl("REF", final_output_wide$snv_read_specific_1),
        paste0("discordant"),
    ifelse(
        grepl("REF", final_output_wide$snv_read_specific_1) &
        grepl("^(?!.*REF).*FALSE.*$",
              final_output_wide$snv_read_specific_2, perl = TRUE),
        paste0("discordant"),
    ifelse(
        grepl("REF", final_output_wide$snv_read_specific_2) &
        grepl("^(?!.*REF).*FALSE.*$",
              final_output_wide$snv_read_specific_1, perl = TRUE),
        paste0("discordant"),
    ifelse(
        grepl("REF", final_output_wide$snv_read_specific_1) &
        grepl("REF", final_output_wide$snv_read_specific_2),
        paste0("REF"),
    ifelse(
        grepl("REF", final_output_wide$snv_read_specific_2) &
        grepl("REF", final_output_wide$snv_read_specific_1),
        paste0("REF"),
    ifelse(
        grepl("^(?!.*REF).*FALSE.*$",
            final_output_wide$snv_read_specific_1, perl = TRUE) &
        is.na(final_output_wide$snv_read_specific_2),
        paste0(gsub(":([^:]+)$", "",
               final_output_wide$snv_read_specific_1), sep = ":not_mut"),
    ifelse(
        grepl("^(?!.*REF).*FALSE.*$",
              final_output_wide$snv_read_specific_2, perl = TRUE) &
        is.na(final_output_wide$snv_read_specific_1),
        paste0(gsub(":([^:]+)$", "", final_output_wide$snv_read_specific_2),
               sep = ":not_mut"),
    ifelse(
        grepl("^(?!.*REF).*FALSE.*$",
            final_output_wide$snv_read_specific_1, perl = TRUE) &
        grepl("^(?!.*REF).*FALSE.*$",
            final_output_wide$snv_read_specific_2, perl = TRUE),
        paste0(gsub(":([^:]+)$", "",
               final_output_wide$snv_read_specific_1), sep = ":not_mut"),
    ifelse(
        grepl("TRUE", final_output_wide$snv_read_specific_1) &
        grepl("TRUE", final_output_wide$snv_read_specific_2),
        gsub(":([^:]+)$", "", final_output_wide$snv_read_specific_1),
  NA))))))))))))))

  ############################################################################
  # Summarise fragment by R1/R2 mismatch and concordance
  ############################################################################
  # Extract the character after the second colon from snv_read_specific_1
  base1 <- sub("^[^:]*:[^:]*:([^:]).*", "\\1", final_output_wide$snv_read_specific_1)

  # Extract the character after the second colon from snv_read_specific_2
  base2 <- sub("^[^:]*:[^:]*:([^:]).*", "\\1", final_output_wide$snv_read_specific_2)

  final_output_wide$mismatch_status  <- ifelse(
        grepl("TRUE", final_output_wide$concordance) &
        !grepl("REF", final_output_wide$match),
        paste0(final_output_wide$match, sep = ":concordant"),
    ifelse(
        is.na(final_output_wide$concordance) &
        !grepl("REF", final_output_wide$match),
        paste0(final_output_wide$match, sep = ":single_read"),
    ifelse(
        is.na(final_output_wide$concordance) &
        grepl("REF", final_output_wide$match),
        paste0(final_output_wide$bed_locus, sep = ":REF:single_read"),
    ifelse(
        !is.na(final_output_wide$concordance) &
        grepl("REF", final_output_wide$match),
        paste0(final_output_wide$bed_locus, sep = ":REF:concordant"),
    paste0(final_output_wide$bed_locus, sep = ":", base1, base2, sep = ":discordant")))))

  ############################################################################
  # If fragment has only one read mutated from the pair, assign strand info
  ############################################################################
  final_output_wide$single_read_strand <- ifelse(
    grepl("single_read", final_output_wide$mismatch_status) &
    !is.na(final_output_wide$snv_read_specific_1),
    str_extract(final_output_wide$snv_read_specific_1,
                "(?<=:)[^:|]*(?=\\|[^|]*$)"),
    ifelse(
    grepl("single_read", final_output_wide$mismatch_status) &
    !is.na(final_output_wide$snv_read_specific_2),
    str_extract(final_output_wide$snv_read_specific_2,
                "(?<=:)[^:|]*(?=\\|[^|]*$)"),
    ""))

  ############################################################################
  # Tidy up Dataframe
  ############################################################################
  final_output_wide$fragment_mismatch_bed <- paste0(
    final_output_wide$mismatch_status,
    final_output_wide$single_read_strand)

  final_df <- final_output_wide  %>% select(fragment_id, fragment_mismatch_bed)

  colnames(final_df) <- c("fragment_id", "locus_info")

  final_df$locus <- ifelse(!grepl("not_mut", final_df$locus_info),
                          gsub("(:[^:]+){2}$", "", final_df$locus_info),
                          ifelse(grepl("not_mut", final_df$locus_info),
                          gsub("(:[^:]+){3}$", "", final_df$locus_info), NA))

  ############################################################################
  # Final tidy up and slight renaming
  ############################################################################
  final_df$locus_status <- ifelse(
        !grepl("REF", final_df$locus_info) &
        !grepl("not_mut", final_df$locus_info),
        paste0("MUT", sep = ":", gsub(".*:([^:]+)$", "\\1",
               final_df$locus_info)),
    ifelse(
        grepl("REF", final_df$locus_info),
        paste0("REF", sep = ":",
               gsub(".*:([^:]+)$", "\\1", final_df$locus_info)),
    ifelse(
        grepl("not_mut", final_df$locus_info),
        paste0("other_base", sep = ":",
               gsub(".*:([^:]+)$", "\\1", final_df$locus_info)),
    paste0("discordant"))))

  ############################################################################
  # Remove redundant strand information
  ############################################################################
  final_df$locus_status <- gsub("[+|-]", "", final_df$locus_status)

  #############################################################################
  # process fragment mismatch dataframe for merging
  #############################################################################
  mismatch_df <- final_df  %>% select(locus_info,
                                      locus_status,
                                      fragment_id)

  rownames(mismatch_df) <- final_df$fragment_id

  return(mismatch_df)
}



#' Write Mutation Table
#'
#' This function processes a GRanges object, summarizing cfDNA fragment
#' information for each target mutation locus.
#' It annotates each locus with the number and type of supporting fragments.
#' Each mismatch type is annotated with the median fragment length.
#' Consensus mismatch is determined by selecting the most frequent mismatch,
#' giving priority to the target mutation's ALT base over other bases.
#' Consensus mismatch is used to derive the trinucleotide substitution (SBS96).
#' The resulting table is written to a .tsv file.
#'
#' @import GenomicRanges
#' @import dplyr
#' @import BSgenome
#' @import IRanges
#' @import stringr
#'
#' @param gr GRanges object containing genomic ranges and associated data.
#' @param output_dir Character string specifying the directory where the output file will be saved.
#' @param output_file Character string specifying the name of the output file (without extension).
#'
#' @return The path to the output file as a character string.
#'
#' @examples
#' output_path <- writeMutTable(gr, "output_directory", "output_filename")
#' print(paste("File written to:", output_path))
writeMutTable <- function(gr, output_dir, output_file) {
  # Ensure the output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  trinuc_df <- callTrinucleotide(gr = gr)

  # Write the dataframe to a .tsv file
  output_path <- file.path(output_dir, paste0(output_file, ".tsv"))
  write.table(trinuc_df, file = output_path, sep = "\t", row.names = FALSE, quote = FALSE)

  return(output_path)
}

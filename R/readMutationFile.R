readMutationFile <- function(mutation_file) {

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
      stringsAsFactors = FALSE, quote = "")
  } else {
    bed <- read.table(mutation_file,
      header = FALSE,
      sep = "\t",
      stringsAsFactors = FALSE, quote = "",
      col.names = c("chr", "pos", "ref", "alt"))
  }

  # Check whether the file format is correct
  result <- check_columns(bed)

  if (result == "Mutation file OK!") {
    message("Mutation file OK!")
    # Process the file
    bed["start"] <- bed["pos"]
    bed["end"] <- bed["pos"]
    bed_full <- bed  %>% select(chr, start, end, ref, alt)

  return(bed_full)

  # Error message
  } else {
      message("Incorrect Mutation File Format")
      message("Please provide a tab separated file with 4 columns:")
      message("chr	pos	ref	alt")
      message("Column header is not required")
  }
}

# Internal Function to check the mutation file
check_columns <- function(df) {

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

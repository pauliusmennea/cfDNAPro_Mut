% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internalMutFunctions.R
\name{readMutations}
\alias{readMutations}
\title{Read and Process the provided Mutation File}
\usage{
readMutations(
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
  frag
)
}
\arguments{
\item{frag}{}
}
\value{
Returns a GRanges object with curated mutational information
}
\description{
Read and Process the provided Mutation File
}

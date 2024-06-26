% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/callTrinucleotide.R
\name{callTrinucleotide}
\alias{callTrinucleotide}
\title{Call trinucleotides and summarise the cfDNA information
for each target mutation locus.}
\usage{
callTrinucleotide(gr)
}
\arguments{
\item{gr}{GRanges object containing genomic ranges and associated data.}
}
\value{
dataframe with summarised mutational and trinucleotide data
}
\description{
This function processes a GRanges object, summarizing cfDNA fragment
information for each target mutation locus.
It annotates each locus with the number and type of supporting fragments.
Each mismatch type is annotated with the median fragment length.
Consensus mismatch is determined by selecting the most frequent mismatch,
giving priority to the target mutation's ALT base over other bases.
Consensus mismatch is used to derive the trinucleotide substitution (SBS96).
}
\examples{
trinuc_df <-  callTrinucleotide(gr)
}

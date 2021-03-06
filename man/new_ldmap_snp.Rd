% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/snp.R
\name{new_ldmap_snp}
\alias{new_ldmap_snp}
\title{create a new ldmap_snp vector}
\usage{
new_ldmap_snp(
  chrom = integer(),
  pos = integer(),
  ref = new_ldmap_allele(),
  alt = new_ldmap_allele()
)
}
\arguments{
\item{chrom}{a vector that can be coerced to integer through
\code{as_chromosome} (e.g c(1L), c(12.0), c("chrX","chr11"))}

\item{pos}{a vector that can be coerced to integer}

\item{ref}{a vector of single nucleotide reference alleles}

\item{alt}{a vector of single nucleotide alternate alleles}
}
\value{

}
\description{
create a new ldmap_snp vector
}

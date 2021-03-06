% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plink.R
\name{bim_cols}
\alias{bim_cols}
\title{Column specification for a plink bim file}
\usage{
bim_cols(
  chrom = col_chromosome(prefix_chr = TRUE),
  rsid = readr::col_character(),
  map = readr::col_double(),
  pos = readr::col_integer(),
  alt = readr::col_character(),
  ref = readr::col_character(),
  ...
)
}
\arguments{
\item{chrom}{chromosome specifier}

\item{rsid}{snp id specifier}

\item{map}{genetic map specifier}

\item{pos}{position specifier}

\item{alt}{alt allele specifier}

\item{ref}{ref allele specifier}

\item{...}{unused}
}
\value{
\code{readr::cols} specification
}
\description{
Column specification for a plink bim file
}

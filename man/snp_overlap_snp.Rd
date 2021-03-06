% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{snp_overlap_snp}
\alias{snp_overlap_snp}
\title{Check SNPs for positional equality}
\usage{
snp_overlap_snp(x, y)
}
\arguments{
\item{x}{vector of query ldmap_snps}

\item{y}{vector of target ldmap_snps}
}
\value{
a vector of integers of length \code{x} that indexes into \code{y}
}
\description{
Check SNPs for positional equality
}

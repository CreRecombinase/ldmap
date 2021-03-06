% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/snp.R
\name{set_chromosomes}
\alias{set_chromosomes}
\title{Change chromosome of a SNP or region}
\usage{
set_chromosomes(x, to)
}
\arguments{
\item{x}{source vector of SNPs}

\item{to}{new chromosome value for SNPs}
}
\value{
a vector of length equal to \code{x} with all(positions(x)==to)
}
\description{
Change chromosome of a SNP or region
}
\examples{
#generate random snps and then set chromosmes to 1:10 and all positions to 10
x <- ldmap_snp(double(10)) \%>\% set_chromosomes(1:10) \%>\% set_positions(10)
all(positions(x)==10)
}

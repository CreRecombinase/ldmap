% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{strand_flip}
\alias{strand_flip}
\title{Determine whether 2 alleles are compatible}
\usage{
strand_flip(ref_alt, reverse = FALSE)
}
\arguments{
\item{query_ref_alt}{is a ref/alt pair}

\item{target_ref_alt}{is another ref/alt pair}
}
\value{
returns a vector with 1 if the query matches the target, -1 if a flip is required, or 0 if they are incompatible;
}
\description{
Determine whether 2 alleles are compatible
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/snp.R
\name{is_ldmap_snp}
\alias{is_ldmap_snp}
\title{Whether or not your vector is a ldmap_snp}
\usage{
is_ldmap_snp(x)
}
\arguments{
\item{x}{vector}
}
\value{
length one boolean
}
\description{
Whether or not your vector is a ldmap_snp
}
\examples{
#FALSE
is_ldmap_snp(letters)
#TRUE
is_ldmap_snp(new_ldmap_snp(1:22,1:100))
}

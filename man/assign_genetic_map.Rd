% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/block.R
\name{assign_genetic_map}
\alias{assign_genetic_map}
\title{Assign or Interpolate Genetic Map Values}
\usage{
assign_genetic_map(snp_df, map_df, strict = FALSE)
}
\arguments{
\item{snp_df}{dataframe of snp coordinates. \code{snp_df} must contain (integer-valued) columns
named \code{chr} and \code{pos}.}

\item{map_df}{dataframe of reference genetic map values. \code{map_df} must contain integer-valued
columns named \code{chr} and \code{pos},
as well as a numeric-valued column called \code{map}, which must be a \emph{strictly} sorted vector of
cumulative genetic map values.}

\item{strict}{boolean indicating whether}
}
\value{
a modified \code{snp_df} with an additional column giving the interpolated genetic map values
at each locus
}
\description{
Assign or Interpolate Genetic Map Values
}

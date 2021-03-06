% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/map.R
\name{jitter_map}
\alias{jitter_map}
\title{Jitter genetic map}
\usage{
jitter_map(map, min_diff = 10 * .Machine$double.eps)
}
\arguments{
\item{map}{cumulative genetic map values.  must by sorted (but need not be strictly sorted)}
}
\value{
numeric vector of length \code{p} with genetic map values
}
\description{
add jitter to cumulative genetic map values so that the results are stricrly sorted
}
\author{
Nicholas Knoblauch
}

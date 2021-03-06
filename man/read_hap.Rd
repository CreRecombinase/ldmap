% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hap.R
\name{read_hap}
\alias{read_hap}
\title{Read \code{.hap} (or \code{.hap.gz}) files into a list of ldmap_ht}
\usage{
read_hap(x, subset = NULL, ...)
}
\arguments{
\item{x}{path to \code{.hap} (or \code{.hap.gz}) file}

\item{subset}{either \code{NULL} ,
in which case all samples are written,
or an integer vector of rows to subset}
}
\value{
list_of ldmap_ht
}
\description{
Read \code{.hap} (or \code{.hap.gz}) files into a list of ldmap_ht
}

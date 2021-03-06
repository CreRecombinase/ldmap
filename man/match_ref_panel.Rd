% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/block.R
\name{match_ref_panel}
\alias{match_ref_panel}
\title{Match Summary Statistics with an external reference panel}
\usage{
match_ref_panel(
  gwas_df,
  ref_snp_struct,
  rsid = integer(),
  remove_ambig = FALSE,
  snp_struct_col = snp_cols(gwas_df),
  flip_sign_col = c("beta"),
  rename_col = TRUE,
  remove_missing = FALSE
)
}
\arguments{
\item{gwas_df}{gwas dataframe}

\item{ref_snp_struct}{ldmap_snp vector of SNP coordinates for the reference panel}

\item{rsid}{optional vector with (integer coded)}

\item{remove_ambig}{whether to remove strand ambiguous query SNPs prior to trying to match target}

\item{snp_struct_col}{column name of ldmap_snp in \code{gwas_df}}

\item{flip_sign_col}{columns that should be flipped when the query SNP is flipped against the target}
}
\value{

}
\description{
Match Summary Statistics with an external reference panel
}

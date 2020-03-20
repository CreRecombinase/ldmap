
#' Formatting method for ldmap allele
#'
#' @param x ldmap_ht
#' @param ... unused
#' @return
#' @export
#' @export format.ldmap_ht
#' @method format ldmap_ht
format.ldmap_ht <- function(x,...){

    ldmap:::format_ht(x)
}

#' stupid/lazy function to guess chromosome from a file path
guess_chr <- function(x){
    file_n <- fs::path_file(x)
    tchr <- stringr::str_replace(file_n,".*(chr[0-9XY]+).*","\\1")
    if (tchr %in% chromosome_levels(as_ucsc = TRUE))
        return(tchr)
    return(NA_character_)
}

##' @title read haplegend
##'
##'
##' @param x path to haplegend file
##' @return tibble with 1 ldmap_snp column
##' @export
read_hap_legend <- function(x, chr = guess_chr(x),read_fun=readr::read_delim) {
    legend_cols <- readr::cols(
        id = readr::col_character(),
        position = readr::col_integer(),
        a0 = readr::col_character(),
        a1 = readr::col_character()
        )
    if(fs::path_ext(x)=="gz")
        stopifnot(fs::path_ext(fs::path_ext_remove(x)) == "legend")
    else
        stopifnot(fs::path_ext(x) == "legend")
    xdf <- read_fun(x,
        delim = " ",
        col_names = names(legend_cols$cols),
        col_types = legend_cols, skip = 1
    ) %>%
        dplyr::transmute(id = id, snp = new_ldmap_snp(chr, position, a0, a1))
}
##' @title read haploid sample info##'
##'
##' @param x `.samples` file
##' @return dataframe with sample info
##' @export
read_hap_samples <- function(x){

    samp_cols <- readr::cols(
            sample = readr::col_character(),
            population = readr::col_character(),
            group = readr::col_character(),
            sex = readr::col_integer()
    )
  ret_df <- readr::read_delim(x,
                              delim = " ",
                              col_names = names(samp_cols$cols),
                              col_types = samp_cols,
                              skip = 1)
}



##' @title read `.hap` (or `.hap.gz`) files into a list of ldmap_ht
##'
##'
##' @param x path to `.hap` (or `.hap.gz`) file
##' @param subset either `NULL` ,in which case all samples are written, or an integer vector of rows to subset
##' @param ...
##' @return list_of ldmap_ht
##' @export
read_hap <- function(x, subset = NULL, ...) {
    if (!is.null(subset)) {
          psimf <- readr::ListCallback$new(function(x, index, subset) parse_hap(x[index %in% subset]), subset = subset)
      } else {
          psimf <- readr::ListCallback$new(function(x, index) parse_hap(x))
      }
    purrr::flatten(readr::read_lines_chunked(x, psimf, ...))
}

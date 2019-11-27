##' column specifier for chromosme
##'
##' @param prefix_chr boolean indicating whether or not there is a `chr` prefix
##' @param ... currently unused
##' @return `col` specification (see `?readr::col`)
##' @export
col_chromosome <- function(prefix_chr =  TRUE, ...) {
    if (prefix_chr)
        return(readr::col_factor(paste0("chr",
                                        c(as.character(1:22),
                                          c("X",
                                            "Y")))))
    return(readr::col_integer())
}

##' Column specification for a plink bim file
##'
##' @param chrom chromosome specifier
##' @param rsid snp id specifier
##' @param map genetic map specifier
##' @param pos position specifier
##' @param alt alt allele specifier
##' @param ref ref allele specifier
##' @param ... unused
##' @return `readr::cols` specification
##' @export
bim_cols <- function(
                     chrom = col_chromosome( prefix_chr = TRUE),
                     rsid =  readr::col_character(),
                     map =  readr::col_double(),
                     pos =  readr::col_integer(),
                     alt =  readr::col_character(),
                     ref =  readr::col_character(), ...) {

    return(readr::cols(chrom = chrom,
                       rsid =  rsid,
                       map =  map,
                       pos =  pos,
                       alt =  alt,
                       ref =  ref))
}



##' Default columns for a `bed` (genomic regions) file
##'
##' @export
bed_region_cols <- function(chrom = col_chromosome(),
                            start =  readr::col_integer(),
                            end =  readr::col_integer(),
                            ...) {
    argl <- rlang::list2(...)
    arg_l <- rlang::list2(
                        chrom = chrom,
                        start = start,
                        end = end,
                        !!!argl)
    return(rlang::exec(readr::cols, !!!arg_l))
}


##'
##' Read a plink bim file
##'
##'
##' @param file path to plink bim file, `file` can be any type you can pass to `readr::read_delim`
##' @param compact a boolean indicating whether to compact chrom,pos,ref and alt as `ldmap_snp` (default to TRUE)
##' @param cols column specification (see `?readr::cols` for more info)
##' @param read_fun function to use to read data(must accept `col_names` and `col_types` arguments)
##' @param ...
##'
##' @return a `tibble` with the contents of the `bim` file
##' @export
read_plink_bim <- function(file, compact = TRUE, cols = bim_cols(),read_fun=readr::read_tsv, ...) {

    ret_df <- read_fun(file, col_names = names(cols$cols),col_types=cols)
    if (compact)
        return(compact_snp_struct(ret_df))
    return(ret_df)
}



##'
##' Read a bed (genomic interval, not binary ped) file
##'
##'
##' @param file path to plink bim file, `file` can be any type you can pass to
##' `readr::read_delim`
##' @param compact a boolean indicating whether to compact chrom,pos,ref and alt
##' as `ldmap_snp` (default to TRUE)
##' @param cols column specification
##' @param read_fun function to use to read data(must accept `col_names` and `col_types` arguments)
##' @param ...
##'
##' @return a `tibble` with the contents of the `bed` file
##' @export
read_bed <- function(file, compact = TRUE, cols = bed_region_cols(),read_fun=readr::read_tsv, ...) {
    bcn <- names(cols$cols)
    ret_df <- read_fun(file, col_names = bcn,
                       col_types = cols,skip = 1L)
    if (compact)
        return(compact_ldmap_range(ret_df, chrom = bcn[1], start = bcn[2], end = bcn[3]))
    return(ret_df)
}

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
        return(compact_ldmap_region(ret_df, chrom = bcn[1], start = bcn[2], end = bcn[3]))
    return(ret_df)
}


#' Read HDF5 file with snp info
#'
#' @param h5file 
#' @param ldmr ldmap_region scalar or vector
#' @param datapath datapath contain snp dataframe
#' @param ... other arguments passed to read_df_h5
#' 
#' Requires `EigenH5`
#'
#' @return a dataframe with only SNPs from inside `ldmr`
#' @export
read_snp_region_h5 <- function(h5file,ldmr,datapath="snp",...){
    if(!requireNamespace("EigenH5", quietly = TRUE)){
        stop("Package \"EigenH5\" needed for this function to work, please install it `(remotes::install_github(\"CreRecombinase/EigenH5\")`)",call.=FALSE)
    }
    argl <- rlang::list2(...)
    if("chrom_offset" %in% EigenH5::ls_h5(h5file)){
        chroms <- unique(chromosomes(ldmr))
        co_df <- read_df_h5(h5file,"chrom_offset",subset=chroms)
        if(nrow(co_df)==1){
            ss=seq(from=co_df$offset,length.out = co_df$datasize)
            snp_df <- rlang::exec(EigenH5::read_df_h5,filename=h5file,datapath = datapath,subset=ss,!!!argl) 
            sc <- snp_df[[snp_cols(snp_df)]]
            return(dplyr::filter(snp_df,sc %overlaps%ldmr))
        }
        else{
            return(purrr::pmap_dfr(co_df,function(offset,datasize,...){
                ss=seq(from=offset,length.out = datasize)
                snp_df <- rlang::exec(EigenH5::read_df_h5,filename=h5file,datapath = datapath,subset=ss,!!!argl)
                sc <- snp_df[[snp_cols(snp_df)]]
                return(dplyr::filter(snp_df,sc %overlaps%ldmr))
            }))
        }
    }
    sinfo <- EigenH5::info_h5(h5file,datapaths = EigenH5::ls_h5(h5file,datapath,full_names = TRUE),attr = TRUE)
    snp_col <- dplyr::filter(sinfo,lengths(attributes)==1) %>% 
        dplyr::filter(purrr::map_lgl(attributes,function(x){
        xc <- x[["class"]]
        if(is.null(xc)){
            return(FALSE)
        }
        return("ldmap_snp"%in% xc)
        })) %>% 
        dplyr::pull(name)
    stopifnot(length(snp_col)==1)
    srv <- EigenH5::read_vector_h5(h5file,snp_col)
    sro <- srv %overlaps% sort(ldmr)
    snp_lgl <- which( sro)
    argl$subset <- snp_lgl
    retdf <-    rlang::exec(EigenH5::read_df_h5,filename=h5file,datapath=datapath,!!!argl)
    return(retdf)

}


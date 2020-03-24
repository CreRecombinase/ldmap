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
                     chrom = col_chromosome(prefix_chr = TRUE),
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




read_fam <- function(x,cols=fam_cols()){
    readr::read_delim(x, delim=" ", col_names=names(cols$cols),col_types=cols)
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

    fc <- scan(as.character(file), what = character(), n = 1, quiet = TRUE)
    chrp <- stringr::str_detect(fc,"^chr")
    levp <- all(stringr::str_detect(cols$cols$chrom$levels,"^chr"))
    if (levp != chrp){
        if(chrp){
            cols$cols$chrom$levels <- paste0("chr",cols$cols$chrom$levels)
        }else{
            cols$cols$chrom$levels <- stringr::str_remove(cols$cols$chrom$levels,"chr")
        }
    }
    ret_df <- read_fun(file, col_names = names(cols$cols),col_types=cols)
    if (compact)
        return(compact_snp_struct(ret_df))
    return(ret_df)
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
fam_cols <- function(
                     fid = readr::col_character(),
                     iid =  readr::col_character(),
                     father_id = readr::col_character(),
                     mother_id = readr::col_character(),
                     sex = readr::col_integer(),
                     phenotype  = readr::col_double() , ...) {

    return(readr::cols(fid = fid,
                       iid =  iid,
                       father_id = father_id,
                       mother_id = mother_id,
                       sex = sex,
                       phenotype  = phenotype))
}


##' Create a vector of plink genotypes from a raw byte vector
##'
##' @param x vector of bytes of type raw
##' @param N number of samples in the data(otherwise assumed it is length(x)*4
##' @return a vector of ldmap_gt
##' @export
new_ldmap_gt <- function(x,N=length(x)*4){
    retv <- new_vctr(x,N=N,class="ldmap_gt")
}


##' Create a vector of plink genotypes from a raw byte vector
##'
##' @param x vector of bytes of type raw
##' @param N number of samples in the data(otherwise assumed it is length(x)*4
##' @return a vector of ldmap_gt
##' @export
new_ldmap_gt <- function(x,N=length(x)*4){
    retv <- new_vctr(x,N=N,class="ldmap_gt")
}




##' Read plink bed file into a (relatively) tidy format
##'
##' @param f path to plink `.bed` file
##' @return
##' @export
read_plink_bed <- function(f,return_fam=FALSE){
    bimf <- fs::path_ext_set(f, "bim")
    bim_df <- read_plink_bim(bimf)
    p  <- nrow(bim_df)
    fam_f <- fs::path_ext_set(f, "fam")
    fam_df <- readr::read_delim(fam_f, delim=" ", col_names=names(fam_cols()$cols),col_types=fam_cols())
    N <- nrow(fam_df)
    ## fc <- file(f, open = "rb", raw = TRUE)
    ## fbytes <- readBin(fc, what = raw(), n = 3)
    ## craw <- as.raw(c(0x6c, 0x1b, 0x01))
    ## stopifnot(all.equal(fbytes, craw))
    bim_df$genotypes  <-  read_plink_bed_l(f,p,N)
    return(bim_df)

}



#' Formatting method for ldmap genotype
#'
#' @param x ldmap_gt
#' @param ... unused
#' @return
#' @export
#' @export format.ldmap_gt
#' @method format ldmap_gt
format.ldmap_gt <- function(x,...){

    ldmap:::format_strings(x)
}




#' @export
`[.ldmap_gt` <-  function(x, i, ...) {
    return(gt_subset(x, i))
}


get_N <- function(x){
    attr(x,"N")
}

##' @title Calculate allele frequency
##'
##'
##' @param x either a vector of type `ldmap_gt`, a list_of of type `ldmap_gt`, or a vector for which
##' @param na.rm whether or not to remove missing data
##' @return
##' @export
calc_AF <- function(x, na.rm = FALSE){

    if (inherits(x,"ldmap_gt"))
        return(gt_af(x,na.rm))
    if (inherits(x,"vctrs_list_of") && inherits(attr(x,"ptype"),"ldmap_gt"))
        return(gt_afs(x,na.rm))
    return(sum(x,na.rm=na.rm) / (2 * length(x[!is.na(x)])))
}

sum_N <- function(x,y){
    stopifnot(inherits(x,"ldmap_gt"),
              inherits(y,"ldmap_gt"))
    get_N(x)+get_N(y)
}

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



#read_plink_fam <- function(


read_plink_bed <- function(f){
    bimf <- fs::path_ext_set(f, "bim")
    bim_df <- read_plink_bim(bimf)
    p  <- nrow(bim_df)
    fam_f <- fs::path_ext_set(f,"fam")
  #  fam_df <-


    if(!inherits(f,"connection")){
        f <- file(f, "rb")
    }
    fbytes <- readBin(f, what = raw(), n = 3)
    craw <- as.raw(c(0x6c, 0x1b, 0x01))
    stopifnot(all.equal(fbytes, craw))
    p <- 141123
    N <-  489
    blocks <- p*ceiling(N/4)
    retvec <- readBin(f,what=raw(),n=blocks)


}

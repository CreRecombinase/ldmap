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



##' @title read plink fam file
##'
##'
##' @param x path to file
##' @param cols fam columns (a readr cols object)
##' @return a dataframe with fam info
##' @export
read_fam <- function(x,cols=fam_cols()){
    readr::read_delim(fs::path_ext_set(x,"fam"), delim=" ", col_names=names(cols$cols),col_types=cols)
}


##' @title read plink fam file
##'
##'
##' @param x path to file
##' @param cols fam columns (a readr cols object)
##' @return a dataframe with fam info
##' @export
read_plink_fam <- function(x, cols = fam_cols()) {
        read_fam(fs::path_ext_set(x,"fam"), cols)
}


##' @title write plink fam data
##'
##'
##' @param fam_df fam dataframe
##' @param file path to write
##' @export
write_plink_fam <- function(fam_df,file){
    readr::write_delim(fam_df,path=fs::path_ext_set(file,"fam"), delim=" ",col_names = FALSE)
}



##' @title write trio of bim,bed and fam files to plink format
##'
##'
##' @param bim_df dataframe with bim info
##' @param bed_l list of ldmap_gt vectors
##' @param fam_df data for individuals
##' @param file file to write (can contain `.bed,.bim or .fam` prefix)
##' @export
write_plink <- function(bim_df, bed_l = bim_df$genotypes, fam_df, file) {
        bimf <- fs::path_ext_set(file, "bim")
        bim_df <- write_plink_bim(bim_df, bimf)
        fam_f <- fs::path_ext_set(file, "fam")
        write_plink_fam(fam_df, fam_f)
        bed_f <- fs::path_ext_set(file, "bed")
        write_plink_bed(bed_l, bed_f)
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

##' @title Write snp dataframe as plink bim file
##'
##'
##' @param bim_df dataframe (as from read_bim) extra columns are ignored
##' @param path output file (or prefix)
##' @param snp_col which column should be expanded (contains ldmap_snp)
##' @param id_col which column has snp names (rsids)
##' @param map_col which column has genetic map values
##' @export
write_plink_bim <- function(bim_df,path,snp_col=snp_cols(bim_df),id_col="rsid",map_col="map",append=FALSE){
    oc <- colnames(bim_df)
    rdf <- explode_snp_struct(bim_df,ldmap_snp=snp_col,alleles_to_character=TRUE,remove=TRUE)
    nc <- colnames(rdf)
    new_c  <- nc[!nc %in% oc]
    alt_col <- new_c[stringr::str_detect(new_c,"alt")]
    ref_col <- new_c[stringr::str_detect(new_c,"ref")]
    stopifnot(!is.null(rdf[[id_col]]))
    rdf[[id_col]] <- int2rsid(rdf[[id_col]])
    if(is.null(rdf[[map_col]]))
        rdf[[map_col]] <- rep(0.0,nrow(rdf))
    readr::write_tsv(rdf[,c("chrom",id_col,map_col,"pos",alt_col,ref_col)],fs::path_ext_set(path,"bim"),col_names = FALSE,append=append)
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

##' @title write listof
##'
##'
##' @param bed_l
##' @param file
##' @return
##' @export
write_plink_bed <- function(bed_l,file,append=FALSE){
    stopifnot(is.list(bed_l))
    ## stopifnot(all(purrr::map_lgl(bed_l,~inherits(.x,"ldmap_gt"))))
    write_plink_bed_f(bed_l, fs::path_ext_set(file,".bed"),append)
    ## out_f <- file(file,open='wb')
    ## writeBin(as.raw(c(0x6c, 0x1b, 0x01)),con=out_f)
    ## purrr::walk(bed_l,~writeBin(vec_data(.x),con=out_f))
}

##' Read plink bed file into a (relatively) tidy format
##'
##' @param f path to plink `.bed` file
##' @return
##' @export
read_plink <- function(f,return_fam=FALSE){
    bimf <- fs::path_ext_set(f, "bim")
    bim_df <- read_plink_bim(bimf)
    p  <- nrow(bim_df)
    fam_f <- fs::path_ext_set(f, "fam")
    fam_df <- readr::read_delim(fam_f, delim=" ", col_names=names(fam_cols()$cols),col_types=fam_cols())
    N <- nrow(fam_df)
    if(return_fam)
        stop("read_plink with return_fam=TRUE is not implemented, just use read_plink with return_fam=FALSE and use `read_plink_fam` to get fam")

    ## fc <- file(f, open = "rb", raw = TRUE)
    ## fbytes <- readBin(fc, what = raw(), n = 3)
    ## craw <- as.raw(c(0x6c, 0x1b, 0x01))
    ## stopifnot(all.equal(fbytes, craw))
    bim_df$genotypes  <-  read_plink_bed_l(f,p,N)
    return(bim_df)
}

##' @title read plink bed file
##'
##'
##' @param f path to bed file (or prefix)
##' @param N number of samples (optional, will read the fam file to figure it out)
##' @param p number of variants (optional will read the bim file to figure it out)
##' @return list_of ldmap_gt
##' @export
read_plink_bed <- function(f,N=NA,p=NA,subset=NULL){
    if(is.na(N))
        N <- nrow(read_plink_fam(fs::path_ext_set(f,"fam")))
    if(!is.null(subset))
        return(read_plink_bed_idx(fs::path_ext_set(f,"bed"),as.integer(subset),N))

    if(is.na(p))
        p <- nrow(read_plink_bim(fs::path_ext_set(f,"bim")))
    read_plink_bed_l(fs::path_ext_set(f,"bed"),p,N)
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

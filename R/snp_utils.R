

#'
#' @method is.na ldmap_snp
#' @export
#' @export is.na.ldmap_snp
is.na.ldmap_snp <- function(x,...){
  is.na(chromosomes(x)) | positions(x)==0
}


#'
#' @method vec_duplicate_id ldmap_snp
#' @export
#' @export vec_duplicate_id.ldmap_snp
#'
vec_duplicate_id.ldmap_snp <- function(x,...){
  vec_duplicate_id(vctrs::vec_data(x))
}



#'
#' @method vec_duplicate_id ldmap_snp
#' @export
#' @export vec_duplicate_id.ldmap_snp
#'
vec_duplicate_id.ldmap_snp <- function(x,...){
  vec_duplicate_id(vctrs::vec_data(x))
}


#'
#' @method vec_duplicate_id ldmap_snp
#' @export
#' @export vec_duplicate_id.ldmap_snp
#'
vec_duplicate_id.ldmap_snp <- function(x,...){
  vec_duplicate_id(vctrs::vec_data(x))
}

#'
#' @method vec_proxy_equal ldmap_snp
#' @export
#' @export vec_proxy_equal.ldmap_snp
vec_proxy_equal.ldmap_snp <- function(x, ...) {
  ldmap_snp_2_dataframe(x)
}



#'
#' @method vec_proxy_equal ldmap_allele
#' @export
#' @export vec_proxy_equal.ldmap_allele
#'
vec_proxy_equal.ldmap_allele <- function(x, ...) {
  as_integer_ldmap_allele(x)
}



#' @method vec_ptype2 ldmap_snp
#' @export
#' @export vec_ptype2.ldmap_snp
#'
vec_ptype2.ldmap_snp <- function(x, y, ...) UseMethod("vec_ptype2.ldmap_snp", y)



#' @method vec_ptype2.ldmap_snp ldmap_snp
#' @export
vec_ptype2.ldmap_snp.ldmap_snp <- function(x, y, ...) new_ldmap_snp()

#' @method vec_ptype2.ldmap_snp double
#' @export
vec_ptype2.ldmap_snp.double <- function(x, y, ...) double()

#' @method vec_ptype2.double ldmap_snp
#' @export
vec_ptype2.double.ldmap_snp <- function(x, y, ...) double()






#' @param from A vector of type ldmap_snp or ldmap_region
#'
#' @param to  A vector of type ldmap_snp or ldmap_region
#' 
#' @return a vector of length equal to the 
#' max(length(from),length(to)) 
#' giving the (signed) distance (in base pairs) from `from` to `to`.
#' Cooridnates on different chromosomes have the maximum distance
#' @export
distance <- function(from, to) {
  UseMethod("distance",from)
}

#' @export distance.ldmap_snp
#' @method distance ldmap_snp
distance.ldmap_snp <- function(from, to) {
  UseMethod("distance.ldmap_snp",to)
}


#' Change position of a SNP
#'
#' @param x source vector of SNPs
#' @param to new position value for SNPs
#'
#' @return a vector of length equal to `x` with all(positions(x)==to)
#' @export
#'
#' @examples
#' x <- ldmap_snp(double(10)) %>% set_chromosomes(1:10) %>% set_positions(10)
#' all(positions(x)==10)
set_positions <- function(x,to){
  y <- x
  positions(y) <- to
  return(y)
}


#' Change chromosome of a SNP or region
#'
#' @param x source vector of SNPs
#' @param to new chromosome value for SNPs
#'
#' @return a vector of length equal to `x` with all(positions(x)==to)
#' @export
#'
#' @examples
#' #generate random snps and then set chromosmes to 1:10 and all positions to 10
#' x <- ldmap_snp(double(10)) %>% set_chromosomes(1:10) %>% set_positions(10)
#' all(positions(x)==10)
set_chromosomes <- function(x,to){
  y <- x
  chromosomes(y) <- to
  return(y)
}





#' @export vec_cast.ldmap_snp
#' @method vec_cast ldmap_snp
#' @export
vec_cast.ldmap_snp <- function(x, to, ...) {
  UseMethod("vec_cast.ldmap_snp")
}


#' @export vec_cast.ldmap_snp.default
#' @method vec_cast.ldmap_snp default
#' @export
vec_cast.ldmap_snp.default <- function(x, to, ...) {
  vec_default_cast(x, to)
}



#' Formatting method for ldmap snps
#'
#' @param x ldmap_snp
#' @param ... unused
#'
#' @return
#' @export
#' @export format.ldmap_snp
#' @method format ldmap_snp
format.ldmap_snp <- function(x, ...) {
  format_ldmap_snp(x)
}




#' @export vec_cast.ldmap_snp.ldmap_snp
#' @export
#' @method vec_cast.ldmap_snp ldmap_snp
vec_cast.ldmap_snp.ldmap_snp <- function(x, to, ..., x_arg = "", to_arg = "") x

#' @export vec_cast.ldmap_snp.double
#' @export
#' @method vec_cast.ldmap_snp double
vec_cast.ldmap_snp.double <- function(x, to, ..., x_arg = "", to_arg = "")
    structure(vctrs::vec_data(x), class = c("ldmap_snp", "vctrs_vctr"))


#' @export vec_cast.double.ldmap_snp
#' @export
#' @method vec_cast.double ldmap_snp
vec_cast.double.ldmap_snp <- function(x, to, ..., x_arg = "", to_arg = "")
    vctrs::vec_data(x)


#' @export vec_cast.integer64.ldmap_snp
#' @export
#' @method vec_cast.integer64 ldmap_snp
vec_cast.integer64.ldmap_snp <- function(x, to, ..., x_arg = "", to_arg = "")
    structure(unclass(x),class="integer64")


#' @export vec_cast.character.ldmap_snp
#' @export
#' @method vec_cast.character ldmap_snp
vec_cast.character.ldmap_snp <- function(x, to, ..., x_arg = "", to_arg = "")
    format_ldmap_snp(x)



#' @export
as_ldmap_snp <- function(x) {
  vec_cast(x, new_ldmap_snp())
}




#' Comparison method
#'
#' @param x a vector of ldmap_snps
#' @param ... unused
#'
#' @return a ranking of ldmap_snps
#'
#' @method vec_proxy_compare ldmap_snp
#' @export
#' @export vec_proxy_compare.ldmap_snp
vec_proxy_compare.ldmap_snp <- function(x, ...) {
  structure(unclass(x),class="integer64")
}


##' Return the column(s) containing ldmap_snps
##'
##' @param df a dataframe
##' @return a character vector with the names `ldmap_snp` columns
##' @export
##' @author Nicholas Knoblauch
snp_cols <- function(df) {
    stopifnot(is.data.frame(df))
    colnames(df)[purrr::map_lgl(df, ~ inherits(.x, "ldmap_snp"))]
}


#' convert snp struct back into seperate columns
#'
#' @param df dataframe with at least one ldmap_snp_struct
#' @param ldmap_snp
#'
#' @return
#' @export
#'
#' @examples
#' #TODO
#'
explode_snp_struct <- function(df,
                               ldmap_snp = snp_cols(df),
                               alleles_to_character = FALSE, remove = TRUE) {

    if (length(ldmap_snp) > 1) {
        warning("Multiple ldmap_snp columns detected,\
 and `ldmap_snp` was not specified, using the first: ", ldmap_snp[1])
        ldmap_snp <- ldmap_snp[1]
    }
    stopifnot(length(ldmap_snp) == 1)

    snpc <- df[[ldmap_snp]]
    if (remove) {
    df <- df[, colnames(df) != ldmap_snp]
  }
  return(dplyr::bind_cols(df, ldmap_snp_2_dataframe(snpc, alleles_to_character)))
}


#' convert SNP columns of a dataframe into a compact ldmap_snp
#'
#'
#' @param df dataframe
#' @param chrom chrom column
#' @param pos position column
#' @param ref reference allele column
#' @param alt alt allele column
#' @param remove whether to remove these columns or keep them after compaction
#' @param snp_struct name of new ldmap_snp column
#'
#' @return copy of original datafame with new column `snp_struct`
#' @export
#'
compact_snp_struct <- function(df,
                               chrom = "chrom",
                               pos = "pos",
                               ref = "ref",
                               alt = "alt",
                               snp_struct = "snp_struct",
                               remove = TRUE) {
  stopifnot(all(c(chrom, pos) %in% colnames(df)))

  use_ref <- !is.na(ref)


  use_alt <- !is.na(alt)
  if (use_alt && !use_ref) {
    warning("'ref' was explicitly skipped by being set to NA\
, but alt was specified or left default, `alt` is being skipped as well")
    use_alt <- FALSE
  }


  if (use_ref && (!c(ref) %in% colnames(df))){
      warn(paste0("ref column: ",
                  ref,
                  " not found in colnames(df), skipping both alleles"))
    use_ref <- FALSE
    use_alt <- FALSE
    nref <- nalt <- NA_character_
  }else{
    nref <- ref
    nalt <- alt
  }

  if (use_alt && (!c(alt) %in% colnames(df))) {
      warn(paste0("alt column: ",
                  alt,
                  " not found in colnames(df), skipping alt allele"))
    use_alt <- FALSE
    nalt <- NA_character_
  }else{
    nalt <- alt
  }
  if(use_ref && use_alt)
    df[[snp_struct]] <- new_ldmap_snp(
      chrom = df[[chrom]],
      pos = df[[pos]],
      ref = df[[nref]],
      alt = df[[nalt]]
    )
  if(use_ref && !use_alt)
    df[[snp_struct]] <- new_ldmap_snp(
      chrom = df[[chrom]],
      pos = df[[pos]],
      ref = df[[nref]],
    )
  if(!use_ref && !use_alt)
    df[[snp_struct]] <- new_ldmap_snp(
      chrom = df[[chrom]],
      pos = df[[pos]],
    )
  if(!use_ref && use_alt)
    stop("alt without ref is not supported")

  if (!remove) {
    return(df)
  }
  ucols <- c(chrom, pos, ref, alt)
  ucols <- ucols[!is.na(ucols)]
  return(dplyr::select(df, -dplyr::one_of(ucols)))
}

#' @method vec_ptype2 ldmap_allele
#' @export
#' @export vec_ptype2.ldmap_allele
#'
vec_ptype2.ldmap_allele <- function(x, y, ...) UseMethod("vec_ptype2.ldmap_allele", y)

#' @method vec_ptype2.ldmap_allele ldmap_allele
#' @export
#' @export vec_ptype2.ldmap_allele.ldmap_allele
vec_ptype2.ldmap_allele.ldmap_allele <- function(x, y, ...)
    new_ldmap_allele()

#' @export vec_ptype2.ldmap_allele.character
#' @method vec_ptype2.ldmap_allele character
#' @export
vec_ptype2.ldmap_allele.character <- function(x, y, ...)
    character()

#' @export vec_ptype2.character.ldmap_allele
#' @method vec_ptype2.character ldmap_allele
#' @export
vec_ptype2.character.ldmap_allele <- function(x, y, ...)
    character()

#' @export vec_ptype2.ldmap_allele.integer
#' @method vec_ptype2.ldmap_allele integer
#' @export
vec_ptype2.ldmap_allele.integer <- function(x, y, ...) integer()

#' @export vec_ptype2.integer.ldmap_allele
#' @method vec_ptype2.integer ldmap_allele
#' @export
vec_ptype2.integer.ldmap_allele <- function(x, y, ...) integer()

#' @export vec_cast.ldmap_allele
#' @method vec_cast ldmap_allele
#' @export
vec_cast.ldmap_allele <- function(x, to, ...) {
  UseMethod("vec_cast.ldmap_allele")
}

#' @export
vec_cast.ldmap_allele.default <- function(x, to, ...) {
  vec_default_cast(x, to)
}


#' @export vec_cast.ldmap_allele.ldmap_allele
#' @export
#' @method vec_cast.ldmap_allele ldmap_allele
vec_cast.ldmap_allele.ldmap_allele <- function(x, to, ..., x_arg = "", to_arg = "") x



#' @export vec_cast.ldmap_allele.character
#' @method vec_cast.ldmap_allele character
#' @export
vec_cast.ldmap_allele.character <- function(x, to, ...) {
  new_ldmap_allele(x)
}


#' @export vec_cast.ldmap_allele.integer
#' @method vec_cast.ldmap_allele integer
#' @export
vec_cast.ldmap_allele.integer <- function(x, to, ...) {
  new_ldmap_allele(x)
}


#' @export vec_cast.character.ldmap_allele
#' @method vec_cast.character ldmap_allele
#' @export
vec_cast.character.ldmap_allele <- function(x, to, ...) {
  format_ldmap_allele(x)
}



#' @export vec_cast.integer.ldmap_allele
#' @method vec_cast.integer ldmap_allele
#' @export
vec_cast.integer.ldmap_allele <- function(x, to, ...) {
  as_integer_ldmap_allele(x)
}

#' Formatting method for ldmap allele
#'
#' @param x ldmap_snp
#' @param ... unused
#' @return
#' @export
#' @export format.ldmap_allele
#' @method format ldmap_allele
format.ldmap_allele <- function(x, ...) {
  format_ldmap_allele(x)
}


#' Whether or not your vector is a ldmap_snp
#'
#' @param x vector
#'
#' @return length one boolean
#' @export
#'
#' @examples
#' #FALSE
#' is_ldmap_snp(letters)
#' #TRUE
#' is_ldmap_snp(new_ldmap_snp(1:22,1:100))
is_ldmap_snp <- function(x) {
  inherits(x, "ldmap_snp")
}



#' New ldmap_snp
#'
#' @param x vector
#'
#' @return ldmap_snp vector
#' @export
#'
ldmap_snp <- function(x=double()){
  as_ldmap_snp(x)
}









#' Change chromosome value for a region or SNP
#'
#' @param x ldmap_region or ldmap_snp
#' @param value new chromosome value
#'
#' @return an ldmap_region or ldmap_snp with the appropriate chromsome
#' @export
#'
`chromosomes<-` <- function(x, value) {
  if(inherits(x,"ldmap_snp")){
    return(new_ldmap_snp(chrom = rep(value,len=length(x)),pos = positions(x),ref = ref_alleles(x),alt_alleles(x))    )
  }
  if(inherits(x,"ldmap_region")){
    return(new_ldmap_region(chrom = rep(value,len=length(x)),start = starts(x),end = ends(x)))
  }
  stop("x must be type ldmap_snp or ldmap_region")
}



#' Check for membership of a SNP in a region/region
#'
#' @param ldmap_snp vector of ldmap_snps
#' @param ldmap_region vector of ldmap_regions
#'
#' @return boolean vector of length(length(ldmap_snp))
#' @export
is_snp_in_region <- function(ldmap_snp,ldmap_region){
  
  stopifnot(inherits(ldmap_snp,"ldmap_snp"),
            inherits(ldmap_region,"ldmap_region"))
  return(!is.na(snp_in_region(ldmap_snp,ldmap_region)))
}


#' Change position for a vector of ldmap_snps
#'
#' @param x ldmap_snp
#' @param value new position
#'
#' @return an ldmap_snp vector with all positions at `value`
#' @export
#'
`positions<-` <- function(x, value) {
  if(inherits(x,"ldmap_snp")){
    return(new_ldmap_snp(chrom = chromosomes(x),pos = rep(value,len=length(x)),ref = ref_alleles(x),alt = alt_alleles(x)))
  }
  stop("x must be type ldmap_snp")
}




#' Change ref_allele for a vector of ldmap_snps
#'
#' @param x ldmap_snp
#' @param value new ref_allele
#'
#' @return an ldmap_snp vector with ref_alleles equal to `value`
#' @export
#'
`ref_alleles<-` <- function(x, value) {
  if(inherits(x,"ldmap_snp")){
    return(new_ldmap_snp(chrom = chromosomes(x),pos = positions(x),ref = rep(value,len=length(x)),alt = alt_alleles(alt)))
  }
  stop("x must be type ldmap_snp")
}


#' Change ref_allele for a vector of ldmap_snps
#'
#' @param x ldmap_snp
#' @param value new ref_allele
#'
#' @return an ldmap_snp vector with ref_alleles equal to `value`
#' @export
#'
`alt_alleles<-` <- function(x, value) {
  if(inherits(x,"ldmap_snp")){
    return(new_ldmap_snp(chrom = chromosomes(x),pos = positions(x),ref = ref_alleles(x),alt = rep(value,len=length(x))))
  }
  stop("x must be type ldmap_snp")
}

##' Clear reference and alternate allele for ldmap_snp
##'
##' @param x a vector of ldmap_snps
##' @return a vector of the length of x with chromosome and position set to `N`
clear_alleles <- function(x){
    stopifnot(inherits(x, "ldmap_snp"))
    y <- x
    ref_alleles(y) <- new_ldmap_allele('N')
    alt_alleles(y) <- ref_alleles(y)
    y
}






#' Overlapping interval
#'
#' @param x vector of ldmap_snp or ldmap_region
#' @param y vector of ldmap_snp or ldmap_region
#'
#' @return vector of ldmap_region with overlap between x and y
#' @export
#'
common_region <- function(x,y){
  dplyr::if_else(chromosomes(x)==chromosomes(y),new_ldmap_region(chromosomes(x),pmax(starts(x),starts(y)),pmin(ends(x),ends(y))),new_ldmap_region(NA_integer_,NA_integer_,NA_integer_))
}


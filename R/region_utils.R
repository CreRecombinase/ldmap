
#' translate a vector of ldmap_regions by a constant factor
#'
#' @param x vector of ldmap_region
#' @param offset amount to move begin and end (recycled as necessary)
#'
#' @return vector of ldmap_regions
#' @export
#'
#' @examples
#' # move all regions over by 100 base pairs
#' ldetect_EUR_100 <- translate_ldmr(ldetect_EUR,100L)
#' #move first region by 100 and second by 0
#' ld2 <- translate_ldmr(ldetect_EUR[1:2],c(100L,0L))
translate_ldmr <- function(x,offset=0L){
  new_ldmap_region(chromosomes(x),starts(x)+offset,ends(x)+offset)
}



#' @export distance.ldmap_region
#' @method distance ldmap_region
distance.ldmap_region <- function(from, to) {
  UseMethod("distance.ldmap_region", to)
}


#' Change end position of a region
#'
#' @param x source vector of region
#' @param to new position value for region
#'
#' @return a vector of length equal to `x` with all(ends(x)==to)
#' @export
#'
#' @examples
#' #generate random snps and then set chromosmes to 1:10 and all positions to 10
#' x <- ldmap_snp(double(10)) %>% set_chromosomes(1:10) %>% set_positions(10)
#' all(positions(x)==10)
set_ends <- function(x,to){
  y <- x
  ends(y) <- to
  return(y)
}


#' Change start position of a region
#'
#' @param x source vector of SNPs
#' @param to new position value for SNPs
#'
#' @return a vector of length equal to `x` with all(positions(x)==to)
#' @export
#'
#' @examples
#' #generate random snps and then set chromosmes to 1:10 and all positions to 10
#' x <- ldmap_snp(double(10)) %>% set_chromosomes(1:10) %>% set_positions(10)
#' all(positions(x)==10)
set_starts <- function(x,to){
  y <- x
  starts(y) <- to
  return(y)
}



##' Return the column(s) containing ldmap_regions
##'
##' @param df a dataframe
##' @return a character vector with the names `ldmap_region` columns
##' @author Nicholas Knoblauch
region_cols <- function(df) {
    stopifnot(is.data.frame(df))
    colnames(df)[purrr::map_lgl(df, ~ inherits(.x, "ldmap_region"))]
}


#' convert ldmap_region back in to individual columns
#'
#' @param df dataframe with at least one ldmap_snp_struct
#' @param ldmap_region ldmap_region
#' @param remove boolean indicating whether `ldmap_region`
#' column should be removed
#'
#' @return dataframe with individual columns
#' @export
#'
#'
explode_ldmap_region <- function(df,
                                ldmap_region = region_cols(df),
                                remove = TRUE) {
  if (length(ldmap_region) > 1) {
      warning("Multiple ldmap_region columns detected,\
 and `ldmap_region` was not specified, using the first: ",
 ldmap_region[1])
      ldmap_region <- ldmap_region[1]
  }
  stopifnot(length(ldmap_region) == 1)
  ldr <- df[[ldmap_region]]
  if (remove) {
      df <- df[, colnames(df) != ldmap_region]
  }
  return(dplyr::bind_cols(df, ldmap_region_2_data_frame(ldr)))
}



#' convert  of a dataframe into a compact ldmap_region
#'
#'
#' @param df dataframe
#' @param chrom chromosome column
#' @param start column name for beginning of region
#' @param end column name for end of region
#' @param ldmap_region name for new ldmap_region column
#' @param remove boolean indicating whether or not to remove aforementioned columns
#' @return copy of original datafame with new column `ldmap_region`
#' @export
#'
compact_ldmap_region <- function(df,
                                chrom = "chrom",
                                start = "start",
                                end = "end",
                                ldmap_region = "ldmap_region",
                                remove = TRUE) {
    stopifnot(all(c(chrom, start, end) %in% colnames(df)))


    df[[ldmap_region]] <- new_ldmap_region(
        chrom = df[[chrom]],
        start = df[[start]],
        end =  df[[end]]
    )
    if (!remove) {
        return(df)
    }
    ucols <- c(chrom, start, end)
    ucols <- ucols[!is.na(ucols)]
    return(dplyr::select(df, -dplyr::one_of(ucols)))
}

##' Predicate testing for whether an object is of a type in `ldmap`
##'
##' @param x object to be tested
##' @return scalar logical
is_ldmap <- function(x){
    inherits_any(x,c("ldmap_region","ldmap_snp","ldmap_allele"))
}



#' Whether or not your vector is a ldmap_region
#'
#' @param x vector
#'
#' @return length one boolean
#' @export
#'
#' @examples
#' #FALSE
#' is_ldmap_region(letters)
#' #TRUE
#' is_ldmap_region(new_ldmap_region(1:22,1:100,200:300))
is_ldmap_region <- function(x) {
  inherits(x, "ldmap_region")
}


#' Width of regions
#'
#' @param x an ldmap_region
#'
#' @return width of the vector
#' @export
widths <- function(x){
  stopifnot(inherits(x,"ldmap_region"))
  return(ends(x)-starts(x))
}




#' Midpoints of regions
#'
#' @param x an ldmap_region
#'
#' @return midpoint of the vector
#' @export
midpoints <- function(x){
  stopifnot(is_ldmap(x))
  return(starts(x)+round((ends(x)-starts(x))/2))
}


#' @export
as_ldmap_region <- function(x) {
  vec_cast(x, new_ldmap_region())
}



#' @export vec_cast.ldmap_region
#' @method vec_cast ldmap_region
#' @export
vec_cast.ldmap_region <- function(x, to, ...) {
  UseMethod("vec_cast.ldmap_region")
}


#' @export vec_cast.ldmap_region.default
#' @method vec_cast.ldmap_region default
#' @export
vec_cast.ldmap_region.default <- function(x, to, ...) {
  vec_default_cast(x, to)
}

#' @export vec_cast.ldmap_region.ldmap_snp
#' @export
#' @method vec_cast.ldmap_region ldmap_snp
vec_cast.ldmap_region.ldmap_snp <- function(x, to, ..., x_arg = "", to_arg = "")
    new_ldmap_region(chrom = chromosomes(x),start = positions(x),end = positions(x)+1)



#' @export vec_cast.double.ldmap_region
#' @export
#' @method vec_cast.double ldmap_region
vec_cast.double.ldmap_region <- function(x, to, ..., x_arg = "", to_arg = "")
    vctrs::vec_data(x)


#' @export vec_cast.integer.ldmap_region
#' @export
#' @method vec_cast.integer ldmap_region
vec_cast.integer.ldmap_region <- function(x, to, ..., x_arg = "", to_arg = "")
    as_integer_ldmap_region(x)


#' @export vec_cast.ldmap_region.ldmap_region
#' @export
#' @method vec_cast.ldmap_region ldmap_region
vec_cast.ldmap_region.ldmap_region <- function(x, to, ..., x_arg = "", to_arg = "")
  x



#' @export vec_cast.character.ldmap_region
#' @export
#' @method vec_cast.character ldmap_region
vec_cast.character.ldmap_region <- function(x, to, ..., x_arg = "", to_arg = "")
    format_ldmap_region(x)



#' @export vec_cast.ldmap_region.character
#' @export
#' @method vec_cast.ldmap_region character
vec_cast.ldmap_region.character <- function(x, to, ..., x_arg = "", to_arg = "")
    parse_ldmap_region(x)


#' @export vec_cast.ldmap_region.double
#' @export
#' @method vec_cast.ldmap_region double
vec_cast.ldmap_region.double <- function(x, to, ..., x_arg = "", to_arg = "")
    structure(vctrs::vec_data(x), class = c("ldmap_region", "vctrs_vctr"))


#' Formatting method for ldmap regions
#'
#' @param x ldmap_region
#' @param ...
#'
#' @return
#' @export
#' @export format.ldmap_region
#' @method format ldmap_region
format.ldmap_region <- function(x, ...) {
  format_ldmap_region(x)
}





#' @method vec_ptype2.ldmap_region default
#' @export
#' @export vec_ptype2.ldmap_region.default
vec_ptype2.ldmap_region.default <- function(x, y, ..., x_arg = "x", y_arg = "y") {
  vec_default_ptype2(x, y, x_arg = x_arg, y_arg = y_arg)
}



#' @method vec_ptype2 ldmap_region
#' @export
#' @export vec_ptype2.ldmap_region
vec_ptype2.ldmap_region <- function(x, y, ...) UseMethod("vec_ptype2.ldmap_region", y)



#' New ldmap_region
#'
#' @param x vector
#'
#' @return ldmap_region
#' @export
#'
#' @examples
#' #under the hood the data is a double
#' as.numeric(ldetect_EUR[1])
#' #which means we can convert it back
#' ldmap_region(as.numeric(ldetect_EUR[1]))
ldmap_region <- function(x=double()){
  as_ldmap_region(x)
}




#' @method vec_ptype2.ldmap_region ldmap_region
#' @export
#' @export vec_ptype2.ldmap_region.ldmap_region
vec_ptype2.ldmap_region.ldmap_region <- function(x, y, ...) new_ldmap_region()

#' @method vec_ptype2.ldmap_region double
#' @export
#' @export vec_ptype2.ldmap_region.double
vec_ptype2.ldmap_region.double <- function(x, y, ...) double()



#' @method vec_ptype2.ldmap_region ldmap_snp
#' @export
#' @export vec_ptype2.ldmap_region.ldmap_snp
vec_ptype2.ldmap_region.ldmap_snp <- function(x, y, ...) ldmap_region()


#' @method vec_ptype2.ldmap_snp ldmap_region
#' @export
#' @export vec_ptype2.ldmap_snp.ldmap_region
vec_ptype2.ldmap_snp.ldmap_region <- function(x, y, ...) ldmap_region()




#' @method vec_ptype2.double ldmap_region
#' @export
#' @export vec_ptype2.double.ldmap_region
vec_ptype2.double.ldmap_region <- function(x, y, ...) double()




#' Change start value for ldmap_region
#'
#' @param x ldmap_region
#' @param value new  start
#'
#' @return an ldmap_region with all starts at `value`
#' @export
#'
`starts<-` <- function(x, value) {
  if(inherits(x,"ldmap_region")){
    return(new_ldmap_region(chrom = chromosomes(x),start = rep(value,len=length(x)),end = ends(x)))
  }
  stop("x must be type ldmap_region")
}




#' Change end value for ldmap_region
#'
#' @param x ldmap_region
#' @param value new end
#'
#' @return an ldmap_region with all starts at `value`
#' @export
#'
`ends<-` <- function(x, value) {
  if(inherits(x,"ldmap_region")){
    return(new_ldmap_region(chrom = chromosomes(x),start = starts(x),end =rep(value,len=length(x))))
  }
  stop("x must be type ldmap_region")
}


#' Check for membership of a region in another region
#'
#' @param query vector of ldmap_regions
#' @param target vector of ldmap_regions
#' @param allow_overlap whether partial overlap counts as membership
#'
#' @return boolean vector of length(query)
#' @export
is_region_in_region <- function(query,target,allow_overlap=FALSE){

  stopifnot(inherits(query,"ldmap_region"),
            inherits(target,"ldmap_region"))
  return(!is.na(region_in_region(query,target,allow_overlap)))
}


#' Check for overlap
#'
#' @param x vector of ldmap_snp or ldmap_region
#' @param y vector of ldmap_snp or ldmap_region
#'
#' @return boolean vector
#' @export
#'
overlaps <- function(x,y){
  dplyr::if_else(chromosomes(x)==chromosomes(y),starts(x)<= ends(y) & starts(y) <= ends(x),FALSE)
}


#' @export
`%overlaps%` <- function(x,y){
    UseMethod("%overlaps%", x)
}

#' @export `%overlaps%.ldmap_region`
#' @method `%overlaps%` ldmap_region
`%overlaps%.ldmap_region` <- function(x,y){
    UseMethod("%overlaps%.ldmap_region",y)
}

#' @export `%overlaps%.ldmap_snp`
#' @method `%overlaps%` ldmap_snp
`%overlaps%.ldmap_snp` <- function(x,y){
    UseMethod("%overlaps%.ldmap_snp",y)
}

#' @export `%overlaps%.ldmap_region.ldmap_region`
#' @method `%overlaps%.ldmap_region` ldmap_region
`%overlaps%.ldmap_region.ldmap_region` <- function(x,y){
    !is.na(region_in_region(x,y,TRUE))
}

#' @export `%overlaps%.ldmap_region.ldmap_snp`
#' @method `%overlaps%.ldmap_region` ldmap_snp
`%overlaps%.ldmap_region.ldmap_snp` <- function(x,y){
    !is.na(region_overlap_snp(x,y))
}


#' @export `%overlaps%.ldmap_snp.ldmap_region`
#' @method `%overlaps%.ldmap_snp` ldmap_region
`%overlaps%.ldmap_snp.ldmap_region` <- function(x,y){
    !is.na(snp_in_region(x,y))
}


#' @export `%overlaps%.ldmap_snp.ldmap_snp`
#' @method `%overlaps%.ldmap_snp` ldmap_snp
`%overlaps%.ldmap_snp.ldmap_snp` <- function(x,y){
    !is.na(snp_overlap_snp(x,y))
}





##' Construct new ldmap_region
##'
##' @param chrom chromosome
##' @param start start
##' @param end end
##' @return
##' @author
##' @export
new_ldmap_region <- function(chrom=integer(), start=integer(), end=integer()){
    stopifnot(is.numeric(start),
              is.numeric(end),
              is.factor(chrom)|| is.numeric(chrom))
  if(is.factor(chrom)){
    stopifnot(all.equal(levels(chrom),c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", 
                                        "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", 
                                        "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", 
                                        "chrX", "chrY")))
  }
  p <- max(length(chrom),length(start),length(end))
  nldmap_region(rep(as.integer(chrom),length.out=p),rep(start,length.out=p),rep(end,length.out=p))
}

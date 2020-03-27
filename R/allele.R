#'
#' @method is.na ldmap_allele
#' @export
#' @export is.na.ldmap_allele
is.na.ldmap_allele <- function(x,...){
    vctrs::vec_data(x)==as.raw(0)
}


#' @method vec_ptype2 ldmap_allele
#' @export
#' @export vec_ptype2.ldmap_allele
#'
vec_ptype2.ldmap_allele <- function(x, y, ...) UseMethod("vec_ptype2.ldmap_allele", y)


#' @method vec_ptype2.ldmap_allele ldmap_allele
#' @export
vec_ptype2.ldmap_allele.ldmap_allele <- function(x, y, ...) new_ldmap_allele()


#' @method vec_ptype2.ldmap_allele character
#' @export
vec_ptype2.ldmap_allele.character <- function(x, y, ...) character()


#' @method vec_ptype2.character ldmap_allele
#' @export
vec_ptype2.character.ldmap_allele <- function(x, y, ...) character()



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


#' @export vec_cast.ldmap_allele.default
#' @method vec_cast.ldmap_allele default
#' @export
vec_cast.ldmap_allele.default <- function(x, to, ...) {
  vec_default_cast(x, to)
}



#' @export vec_cast.ldmap_allele.ldmap_allele
#' @export
#' @method vec_cast.ldmap_allele ldmap_allele
vec_cast.ldmap_allele.ldmap_allele <- function(x, to, ..., x_arg = "", to_arg = "") x

#' @export vec_cast.ldmap_allele.raw
#' @export
#' @method vec_cast.ldmap_allele raw
vec_cast.ldmap_allele.raw <- function(x, to, ..., x_arg = "", to_arg = "")
    structure(vctrs::vec_data(x), class = c("ldmap_allele", "vctrs_vctr"))


#' @export vec_cast.raw.ldmap_allele
#' @export
#' @method vec_cast.raw ldmap_allele
vec_cast.raw.ldmap_allele <- function(x, to, ..., x_arg = "", to_arg = "")
    vctrs::vec_data(x)

#' @export vec_cast.character.ldmap_allele
#' @export
#' @method vec_cast.character ldmap_allele
vec_cast.character.ldmap_allele <- function(x, to, ..., x_arg = "", to_arg = "")
    format_ldmap_allele(x)



#' @export vec_cast.ldmap_allele.character
#' @method vec_cast.ldmap_allele character
#' @export
vec_cast.ldmap_allele.character <- function(x, to, ...) {
  new_ldmap_allele(x)
}



#' @export vec_cast.integer.ldmap_allele
#' @method vec_cast.integer ldmap_allele
#' @export
vec_cast.integer.ldmap_allele <- function(x, to, ...) {
  as_integer_ldmap_allele(x)
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



#' @export
as_ldmap_allele <- function(x) {
  vec_cast(x, new_ldmap_allele())
}



#' Formatting method for ldmap allele
#'
#' @param x ldmap_allele
#' @param ... unused
#' @return
#' @export
#' @export format.ldmap_allele
#' @method format ldmap_allele
format.ldmap_allele <- function(x, ...) {
  format_ldmap_allele(x)
}


#' New ldmap_allele
#'
#' @param x vector
#'
#' @return ldmap_allele vector
#' @export
#'
ldmap_allele <- function(x=raw()){
  as_ldmap_allele(x)
}



##' @title Are ref and alt reversed?
##'
##' @param x a factor as returned from `allele_match`
##' @return a logical vector of length `length(x)`
##' @export
is_reversed <- function(x){
    stopifnot(all.equal(levels(x),c("perfect_match",
                                    "reverse_match",
                                    "complement_match",
                                    "reverse_complement_match",
                                    "ambig_match",
                                    "reverse_ambig_match",
                                    "complement_ambig_match",
                                    "reverse_complement_ambig_match")))
    x %in% c("reverse_match",
             "reverse_ambig_match")
}


#'
#' @method vec_proxy_equal ldmap_snp
#' @export
#' @export vec_proxy_equal.ldmap_snp
#'
vec_proxy_equal.ldmap_snp <- function(x, ...) {
  ldmap_snp_2_dataframe(x)
}

#' @method vec_ptype2 ldmap_snp
#' @export
#' @export vec_ptype2.ldmap_snp
#'
vec_ptype2.ldmap_snp <- function(x,y,...)UseMethod("vec_ptype2.ldmap_snp",y)


vec_ptype2.vctrs_percent.default <- function(x, y, ..., x_arg = "x", y_arg = "y") {
  vec_default_ptype2(x, y, x_arg = x_arg, y_arg = y_arg)
}


#' @method vec_ptype2.ldmap_snp ldmap_snp
#' @export
vec_ptype2.ldmap_snp.ldmap_snp <- function(x, y, ...) new_ldmap_snp()

#' @method vec_ptype2.ldmap_snp double
#' @export
vec_ptype2.ldmap_snp.double <- function(x, y, ...) double()

#' @method vec_ptype2.double ldmap_snp
#' @export
vec_ptype2.double.ldmap_snp <- function(x, y, ...) double()



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
  vec_default_cast(x,to)
}


#' @export
#' @method vec_cast.ldmap_snp ldmap_snp
vec_cast.ldmap_snp.ldmap_snp <- function(x, to, ..., x_arg = "", to_arg = "") x

#' @export vec_cast.ldmap_snp.double
#' @export
#' @method vec_cast.ldmap_snp double
vec_cast.ldmap_snp.double <- function(x, to, ..., x_arg = "", to_arg = "") structure(vctrs::vec_data(x),class=c("ldmap_snp","vctrs_vctr"))

#' @export vec_cast.double.ldmap_snp
#' @export
#' @method vec_cast.double ldmap_snp
vec_cast.double.ldmap_snp <- function(x, to, ..., x_arg = "", to_arg = "") vctrs::vec_data(x)



#' @export
as_ldmap_snp <- function(x) {
  vec_cast(x, new_ldmap_snp())
}


#' Comparison method
#'
#' @param x
#' @param ...
#'
#' @return
#'
#' @method vec_proxy_compare ldmap_snp
#' @export
#' @export vec_proxy_compare.ldmap_snp
vec_proxy_compare.ldmap_snp <- function(x, ...) {
  rank.ldmap_snp(x)
}



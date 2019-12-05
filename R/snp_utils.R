
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
#'
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


## make_annots <- function(ldmap_snp,...){
##   other_regions <- rlang::list2(...)
## }

## make_annot <- function(points,regions){
##   return(as.integer(!is.na(snp_in_range(points,regions))))
## }




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






#' Match Summary Statistics with an external reference panel
#'
#' @param gwas_df gwas dataframe
#' @param ref_snp_struct ldmap_snp vector of SNP coordinates for the reference panel
#' @param rsid optional vector with (integer coded)
#' @param remove_ambig whether to remove strand ambiguous query SNPs prior to trying to match target
#' @param snp_struct_col column name of ldmap_snp in `gwas_df`
#' @param flip_sign_col columns that should be flipped when the query SNP is flipped against the target
#'
#' @return
#' @export
#'
match_ref_panel <- function(gwas_df,
                            ref_snp_struct,
                            rsid = integer(),
                            remove_ambig = FALSE,
                            snp_struct_col = snp_cols(gwas_df),
                            flip_sign_col = c("beta")) {
  snp_struct <- gwas_df[[snp_struct_col]]
  stopifnot(!is.null(snp_struct), inherits(snp_struct, "ldmap_snp"))
  if (remove_ambig) {
    gwas_df <- dplyr::filter(gwas_df, !is_strand_ambiguous(snp_struct))
    snp_struct <- gwas_df[[snp_struct_col]]
  }

  gwas_df <- dplyr::arrange(gwas_df, rank.ldmap_snp(snp_struct))
  snp_struct <- gwas_df[[snp_struct_col]]
  match_df <- dplyr::bind_cols(gwas_df, join_snp(snp_struct, ref_snp_struct, rsid = rsid))
  match_df <- dplyr::mutate_at(
                         match_df,
                         dplyr::vars(flip_sign_col),
                         ~ dplyr::if_else(
                                      match_type %in% c("reverse_match", "reverse_ambig_match"),
                                      ., -.))
  return(match_df)
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


#' @export vec_cast.integer.ldmap_snp
#' @export
#' @method vec_cast.integer ldmap_snp
vec_cast.integer.ldmap_snp <- function(x, to, ..., x_arg = "", to_arg = "")
    as_integer_ldmap_snp(x)


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
  rank.ldmap_snp(x)
}




##' Return the column(s) containing ldmap_ranges
##'
##' @param df a dataframe
##' @return a character vector with the names `ldmap_range` columns
##' @author Nicholas Knoblauch
range_cols <- function(df) {
    stopifnot(is.data.frame(df))
    colnames(df)[purrr::map_lgl(df, ~ inherits(.x, "ldmap_range"))]
}



##' Return the column(s) containing ldmap_snps
##'
##' @param df a dataframe
##' @return a character vector with the names `ldmap_snp` columns
##' @author Nicholas Knoblauch
snp_cols <- function(df) {
    stopifnot(is.data.frame(df))
    colnames(df)[purrr::map_lgl(df, ~ inherits(.x, "ldmap_snp"))]
}



#' convert ldmap_range back in to individual columns
#'
#' @param df dataframe with at least one ldmap_snp_struct
#' @param ldmap_range ldmap_range
#' @param remove boolean indicating whether `ldmap_range`
#' column should be removed
#'
#' @return dataframe with individual columns
#' @export
#'
#'
explode_ldmap_range <- function(df,
                                ldmap_range = range_cols(df),
                                remove = TRUE) {
  if (length(ldmap_range) > 1) {
      warning("Multiple ldmap_range columns detected,\
 and `ldmap_range` was not specified, using the first: ",
 ldmap_range[1])
      ldmap_range <- ldmap_range[1]
  }
  stopifnot(length(ldmap_range) == 1)
  ldr <- df[[ldmap_range]]
  if (remove) {
      df <- df[, colnames(df) != ldmap_range]
  }
  return(dplyr::bind_cols(df, ldmap_range_2_data_frame(ldr)))
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
#' TODO
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
  df[[snp_struct]] <- new_ldmap_snp(
    chrom = df[[chrom]],
    pos = df[[pos]],
    ref = df[[nref]],
    alt = df[[nalt]]
  )
  if (!remove) {
    return(df)
  }
  ucols <- c(chrom, pos, ref, alt)
  ucols <- ucols[!is.na(ucols)]
  return(dplyr::select(df, -dplyr::one_of(ucols)))
}


#' convert  of a dataframe into a compact ldmap_range
#'
#'
#' @param df dataframe
#' @param chrom chromosome column
#' @param start column name for beginning of range
#' @param end column name for end of range
#' @param ldmap_range name for new ldmap_range column
#' @param remove boolean indicating whether or not to remove aforementioned columns
#' @return copy of original datafame with new column `ldmap_range`
#' @export
#'
compact_ldmap_range <- function(df,
                                chrom = "chrom",
                                start = "start",
                                end = "end",
                                ldmap_range = "ldmap_range",
                                remove = TRUE) {
    stopifnot(all(c(chrom, start, end) %in% colnames(df)))


    df[[ldmap_range]] <- new_ldmap_range(
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



#' Whether or not your vector is a ldmap_range
#'
#' @param x vector
#'
#' @return length one boolean
#' @export
#'
#' @examples
#' #FALSE
#' is_ldmap_range(letters)
#' #TRUE
#' is_ldmap_range(new_ldmap_range(1:22,1:100,200:300))
is_ldmap_range <- function(x) {
  inherits(x, "ldmap_range")
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
is_ldmap_range <- function(x) {
  inherits(x, "ldmap_snp")
}


#' @export
as_ldmap_range <- function(x) {
  vec_cast(x, new_ldmap_range())
}



#' @export vec_cast.ldmap_range
#' @method vec_cast ldmap_range
#' @export
vec_cast.ldmap_range <- function(x, to, ...) {
  UseMethod("vec_cast.ldmap_range")
}


#' @export vec_cast.ldmap_range.default
#' @method vec_cast.ldmap_range default
#' @export
vec_cast.ldmap_range.default <- function(x, to, ...) {
  vec_default_cast(x, to)
}



#' @export vec_cast.double.ldmap_range
#' @export
#' @method vec_cast.double ldmap_range
vec_cast.double.ldmap_range <- function(x, to, ..., x_arg = "", to_arg = "")
    vctrs::vec_data(x)


#' @export vec_cast.integer.ldmap_range
#' @export
#' @method vec_cast.integer ldmap_range
vec_cast.integer.ldmap_range <- function(x, to, ..., x_arg = "", to_arg = "")
    as_integer_ldmap_range(x)


#' @export vec_cast.ldmap_range.ldmap_range
#' @export
#' @method vec_cast.ldmap_range ldmap_range
vec_cast.ldmap_range.ldmap_range <- function(x, to, ..., x_arg = "", to_arg = "")
  x




#' @export vec_cast.character.ldmap_range
#' @export
#' @method vec_cast.character ldmap_range
vec_cast.character.ldmap_range <- function(x, to, ..., x_arg = "", to_arg = "")
    format_ldmap_range(x)


#' @export vec_cast.ldmap_range.character
#' @export
#' @method vec_cast.ldmap_range character
vec_cast.ldmap_range.character <- function(x, to, ..., x_arg = "", to_arg = "")
    parse_ldmap_range(x)


#' @export vec_cast.ldmap_range.double
#' @export
#' @method vec_cast.ldmap_range double
vec_cast.ldmap_range.double <- function(x, to, ..., x_arg = "", to_arg = "")
    structure(vctrs::vec_data(x), class = c("ldmap_range", "vctrs_vctr"))



#' Formatting method for ldmap ranges
#'
#' @param x ldmap_range
#' @param ...
#'
#' @return
#' @export
#' @export format.ldmap_range
#' @method format ldmap_range
format.ldmap_range <- function(x, ...) {
  format_ldmap_range(x)
}





#' @method vec_ptype2.ldmap_range default
#' @export
#' @export vec_ptype2.ldmap_range.default
vec_ptype2.ldmap_range.default <- function(x, y, ..., x_arg = "x", y_arg = "y") {
  vec_default_ptype2(x, y, x_arg = x_arg, y_arg = y_arg)
}



#' @method vec_ptype2 ldmap_range
#' @export
#' @export vec_ptype2.ldmap_range
vec_ptype2.ldmap_range <- function(x, y, ...) UseMethod("vec_ptype2.ldmap_range", y)


#' New ldmap_range
#'
#' @param x vector
#'
#' @return ldmap_range
#' @export
#'
#' @examples
#' #under the hood the data is a double
#' as.numeric(ldetect_EUR[1])
#' #which means we can convert it back
#' ldmap_range(as.numeric(ldetect_EUR[1]))
ldmap_range <- function(x=double()){
  as_ldmap_range(x)
}

#' @method vec_ptype2.ldmap_range ldmap_range
#' @export 
#' @export vec_ptype2.ldmap_range.ldmap_range
vec_ptype2.ldmap_range.ldmap_range <- function(x, y, ...) new_ldmap_range()

#' @method vec_ptype2.ldmap_range double
#' @export
#' @export vec_ptype2.ldmap_range.double
vec_ptype2.ldmap_range.double <- function(x, y, ...) double()

#' @method vec_ptype2.double ldmap_range
#' @export
#' @export vec_ptype2.double.ldmap_range
vec_ptype2.double.ldmap_range <- function(x, y, ...) double()

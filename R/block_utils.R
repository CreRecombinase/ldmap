#' Assign or Interpolate Genetic Map Values
#'
#' @param snp_df dataframe of snp coordinates. `snp_df` must contain (integer-valued) columns
#' named `chr` and `pos`.
#' @param map_df dataframe of reference genetic map values. `map_df` must contain integer-valued
#' columns named `chr` and `pos`,
#' as well as a numeric-valued column called `map`, which must be a _strictly_ sorted vector of
#' cumulative genetic map values.
#' @param strict boolean indicating whether
#' @return a modified `snp_df` with an additional column giving the interpolated genetic map values
#' at each locus
#' @export
assign_genetic_map <- function(snp_df, map_df, strict = FALSE) {
    ## snp_dfl <- split(snp_df, chromosomes(snp_df$snp))
    ## map_dfl <- dplyr::filter(map_df, u_chr, by = "chr") %>% split(.$chr)
    ## stopifnot(all(names(map_dfl) == names(snp_dfl)))

    ## retdf <- purrr::map2_df(map_dfl,
    ##                         snp_dfl,
    ##                         ~ dplyr::mutate(.y,
    ##                                         map = interpolate_genetic_map(.x$map,
    ##                                                                       positions(.x$snp),
    ##                                                                       positions(.y$snp),
    ##                                                                       strict = strict)))

    return(dplyr::mutate(snp_df,
                         map=new_interpolate_genetic_map(map_df$map,map_df$snp,snp)))
}


calc_theta <- function(m) {
    nmsum <- sum(1 / (1:(2 * m - 1)))
    (1 / nmsum) / (2 * m + 1 / nmsum)
}






#' Query a a reference panel for a set of SNPs and see if any of them need a sign flip
#'
#' @param query_chr chromosome for query SNPs coded as integer or character (coding must match `panel_chr`)
#' @param query_pos position of query SNPs
#' @param query_allele allele for query SNPs
#' @param query_block optional vector specifiyng a chunking scheme (default is NULL)
#' @param panel_chr chromosome for the reference SNPs coded as integer or character (coding must match `query_chr`)
#' @param panel_pos position of reference SNPs
#' @param panel_allele allele for reference SNPs
#' @param panel_block optional vector specifying a chunking scheme (default is NULL)
#'
#' @return an `integer` vector of length equal to `length(query_chr)` where a `1` indicates that the query position
#' is in the reference panel, a `0` indicates that it is not, and a `-1` indicates a match
#' , but that a flip is necessary
#' @export
#'
#' @examples
#'
#' # read in the example data
#' example_bim <- fs::path_package(package = "ldshrink", "test_data/reference_genotype.bim")
#' panel_df <- read.table(example_bim, sep = "	", col.names = c("chr", "snp", "map", "pos", "ref", "alt"))
#' # create "allele" column
#' panel_allele <- paste0(panel_df$ref, ",", panel_df$alt)
#' stopifnot(all(nchar(panel_allele) == 3))
query_reference_snpset <- function(query_chr, query_pos, query_allele, query_block = NULL, panel_chr, panel_pos, panel_allele, panel_block = NULL) {
    stopifnot(is.null(query_block) == is.null(panel_block))

    if ((is.null(panel_block) & is.null(query_block))) {
        join_cols <- c("chr", "pos")
    } else {
        join_cols <- c("chr", "pos", "block")
    }
    if (!is.null(query_block)) {
        query_df <- tibble::tibble(chr = query_chr, pos = query_pos, allele = query_allele, block = query_block)
        panel_df <- tibble::tibble(chr = panel_chr, pos = panel_pos, allele = panel_allele, block = panel_block)
    } else {
        query_df <- tibble::tibble(chr = query_chr, pos = query_pos, allele = query_allele)
        panel_df <- tibble::tibble(chr = panel_chr, pos = panel_pos, allele = panel_allele)
    }
    ij_df <- dplyr::inner_join(query_df, panel_df, by = join_cols)
    ij_df <- dplyr::mutate(ij_df, allele_match = flip_alleles(allele.x, allele.y))
    ij_df <- dplyr::right_join(ij_df, query_df)
    ij_df <- dplyr::mutate(ij_df, allele_match = dplyr::if_else(is.na(allele_match), 0L, allele_match))

    dplyr::pull(ij_df, allele_match)
}



#' Join dataframes with ldmap_snps
#'
#' @param x dataframe with at least one ldmap_snp column to join on
#' @param y dataframe with at least one ldmap_snp column
#' @param by column name of snp column
#' @param suffix suffix to add
#' @param snp_struct_col column name of ldmap_snp in `gwas_df`
#' @param flip_sign_col columns that should be flipped when the query SNP is flipped against the target
#' @param rename_col whether to rename column
#' @param remove_missing whether to remove missing data
#' @param remove_ambig whether to remove strand ambiguous query SNPs prior to trying to match target
#' @return
#' @export
#'
left_snp_join <- function(x,
                          y,
                          by = NULL,
                          suffix = c(".x", ".y"),
                          snp_struct_col = snp_cols(gwas_df),
                          flip_sign_col = c("beta"),
                          rename_col = TRUE,
                          remove_missing = FALSE) {

    if (is.null(by)) {
        snpc_a <- snp_cols(x)
        snpc_b <- snp_cols(b)
    } else {
        stopifnot(length(by) != 1)
        if (is.null(names(by))) {
            snpc_a <- by
        } else {
            snpc_a <- names(by)
        }
        snpc_b <- by
    }
    stopifnot(
        snpc_a %in% colnames(x),
        snpc_b %in% colnames(x)
    )

    asnpcol <- x[[snpc_a]]
    bsnpcol <- x[[snpc_b]]
    stopifnot(inherits(asnpcol,"ldmap_snp"),
              inherits(bsnpcol,"ldmap_snp"))
    jsnp_df <- join_snp(asnpcol,bsnpcol)



}



#' Match Summary Statistics with an external reference panel
#'
#' @param gwas_df gwas dataframe
#' @param ref_snp_struct ldmap_snp vector of SNP coordinates for the reference panel
#' @param rsid optional vector with (integer coded)
#' @param remove_ambig whether to remove strand ambiguous query SNPs prior to trying to match target
#' @param snp_struct_col column name of ldmap_snp in `gwas_df`
#' @param flip_sign_col columns that should be flipped when the query SNP is flipped against the target
#' @return
#' @export
#'
match_ref_panel <- function(gwas_df,
                            ref_snp_struct,
                            rsid = integer(),
                            remove_ambig = FALSE,
                            snp_struct_col = snp_cols(gwas_df),
                            flip_sign_col = c("beta"),
                            rename_col = TRUE,
                            remove_missing = FALSE) {
    snp_struct <- gwas_df[[snp_struct_col]]
    stopifnot(!is.null(snp_struct), inherits(snp_struct, "ldmap_snp"))
    if (remove_ambig) {
        gwas_df <- dplyr::filter(gwas_df, !is_strand_ambiguous(snp_struct))
        snp_struct <- gwas_df[[snp_struct_col]]
    }
    gwas_df <- dplyr::arrange(gwas_df, rank.ldmap_snp(snp_struct))
    snp_struct <- gwas_df[[snp_struct_col]]
    match_df <- dplyr::bind_cols(gwas_df, join_snp(snp_struct, ref_snp_struct, rsid = rsid))
    if (!is.na(flip_sign_col) && is.character(flip_sign_col)) {
        stopifnot(flip_sign_col %in% colnames(gwas_df))
        match_df <- dplyr::mutate_at(
            match_df,
            dplyr::vars(flip_sign_col),
            ~ dplyr::if_else(
                match_type %in% c("reverse_match", "reverse_ambig_match"),
                ., -.
            )
        )
    }
    if (any(is.na(match_df$index))) {
        if (!remove_missing) {
            warning("Some SNPs in gwas_df do not have matching SNPs in ref_snp_struct",
                    " use `remove_missing=TRUE` to drop")
        } else {
            match_df <- dplyr::filter(match_df, !is.na(index))
        }
    }
    if (rename_col) {
        match_df <- dplyr::select(match_df, -{{ snp_struct_col }}) %>%
            dplyr::rename({{ snp_struct_col }} := match)
    }
    return(match_df)
}





#' Estimate LD scores from the LD matrix
#'
#' @param R a pxp LD matrix (as obtained from `estimate_LD`)
#' @param N a scalar representing the number of samples in the reference panel
#'
#' @return A length p vector with LD scores at each locus
#' @export
#'
#' @examples
#' data(reference_genotype)
#' data(reference_map)
#' R <- estimate_LD(reference_panel = reference_genotype, map = reference_map)
#' L2 <- estimate_LDscores(R, nrow(reference_genotype))
estimate_LDscores <- function(R, N) {
    denom <- ifelse(N > 2, N - 2, N)
    if (inherits(R, "Matrix") || inherits(R, "matrix")) {
        L2 <- apply(R^2 - (1 - R^2) / denom, 2, sum)
    } else {
        if (inherits(R, "data.frame")) {
            L2 <- purrr::map_dbl(R$data, ~ sum(.x$r^2 - (1 - .x$r^2) / denom)) - 1
        } else {
            stop("R must be a data.frame,matrix or Matrix")
        }
    }
    return(L2)
}

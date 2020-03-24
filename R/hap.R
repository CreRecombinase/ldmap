
#' Formatting method for ldmap allele
#'
#' @param x ldmap_ht
#' @param ... unused
#' @return
#' @export
#' @export format.ldmap_ht
#' @method format ldmap_ht
format.ldmap_ht <- function(x,...){

    ldmap:::format_ht(x)
}



#' stupid/lazy function to guess chromosome from a file path
guess_chr <- function(x){
    file_n <- fs::path_file(x)
    tchr <- stringr::str_replace(file_n,".*(chr[0-9XY]+).*","\\1")
    if (tchr %in% chromosome_levels(as_ucsc = TRUE))
        return(tchr)
    return(NA_character_)
}

##' @title read haplegend
##'
##'
##' @param x path to haplegend file
##' @return tibble with 1 ldmap_snp column
##' @export
read_hap_legend <- function(x, chr = guess_chr(x),read_fun=readr::read_delim) {
    legend_cols <- readr::cols(
        id = readr::col_character(),
        position = readr::col_integer(),
        a0 = readr::col_character(),
        a1 = readr::col_character()
        )
    if (fs::path_ext(x) == "gz")
        stopifnot(fs::path_ext(fs::path_ext_remove(x)) == "legend")
    else
        stopifnot(fs::path_ext(x) == "legend")
    xdf <- read_fun(x,
        delim = " ",
        col_names = names(legend_cols$cols),
        col_types = legend_cols, skip = 1
    ) %>%
        dplyr::transmute(id = id, snp = new_ldmap_snp(chr, position, a0, a1))
}


##' @title read haploid sample info##'
##'
##' @param x `.samples` file
##' @return dataframe with sample info
##' @export
read_hap_samples <- function(x) {

    samp_cols <- readr::cols(
            sample = readr::col_character(),
            population = readr::col_character(),
            group = readr::col_character(),
            sex = readr::col_integer()
    )
  ret_df <- readr::read_delim(x,
                              delim = " ",
                              col_names = names(samp_cols$cols),
                              col_types = samp_cols,
                              skip = 1)
}


##' @title Read `.hap` (or `.hap.gz`) files into a list of ldmap_ht
##'
##'
##' @param x path to `.hap` (or `.hap.gz`) file
##' @param subset either `NULL` ,
##' in which case all samples are written,
##' or an integer vector of rows to subset
##' @param ...
##' @return list_of ldmap_ht
##' @export
read_hap <- function(x, subset = NULL, ...) {
    subset_factory <- function(subset){
        if (!is.null(subset)) {
            subset <- as.integer(subset)
            stopifnot(all(!is.na(subset)))
            return(function(x, index) {
                idx <- seq(from = index,
                           length.out = length(x))
                parse_hap(x[idx %in% subset])
            })
        } else {
            return(function(x, index) parse_hap(x))
        }
    }
    purrr::flatten(readr::read_lines_chunked(x,readr::ListCallback$new(subset_factory(subset), ...)))
}


##' @title read regions of
##'
##'
##' @param hap_dir directory with haplotype data
##' @param ldmr ldmap_region(s) to read
##' @param hap_file_pattern `glue`-compatible file pattern for haplotype files,
##' when '.legend.gz' is added as a suffix, and `chrom` is substituted, should resolve to a
##' haplotype legend file.  Similarly when `.hap.gz` is added as a suffix should resolve to a
##' haplotype file
##' @return a dataframe or list of dataframes (see `bind_rows`)
##' @export
read_hap_snps <- function(hap_dir,
                            rsid,
                            hap_file_pattern = "EUR.chr{chrom}_1kg_geno") {

    all_files <- fs::dir_ls(hap_dir)
    all_f <- fs::path_file(all_files)
    chromosomes <- 1:22
    legend_cols <- readr::cols(
                              id = readr::col_character(),
                              position = readr::col_integer(),
                              a0 = readr::col_character(),
                              a1 = readr::col_character()
                          )
    callback_factory <- function(chr, snps) {
        function(x, pos) {
            dplyr::transmute(x,
                             id = fast_str2int(id,prefix = "rs"),
                             snp = new_ldmap_snp(chr, position, a0, a1),
                             snp_id = seq(from = pos,
                                          length.out = dplyr::n())) %>%
                dplyr::filter(id %in% snps) %>%
                dplyr::mutate(snp = as.character(snp))
        }
    }
    purrr::map_dfr(chromosomes,
            function(chrom) {
                leg_file <- fs::path(hap_dir,
                                     glue::glue(hap_file_pattern),
                                     ext = "legend.gz")
                hap_file <- fs::path(hap_dir,
                                     glue::glue(hap_file_pattern),
                                     ext = "hap.gz")
                stopifnot(all(file.exists(leg_file)),
                          all(file.exists(hap_file)))
                cbf <- readr::DataFrameCallback$new(callback_factory(chrom, rsid))
                rdf <- readr::read_delim_chunked(file = leg_file,
                                                 delim = " ",
                                                 callback = cbf,
                                                 col_names = names(legend_cols$cols),
                                                 col_types = legend_cols,
                                                 skip = 1)
                dplyr::mutate(rdf, snp = as_ldmap_snp(snp),
                              hap = read_hap(hap_file,
                                             subset = snp_id)) %>%
                    dplyr::select(-snp_id)
            })
}






##' @title read regions of
##'
##'
##' @param hap_dir directory with haplotype data
##' @param ldmr ldmap_region(s) to read
##' @param hap_file_pattern `glue`-compatible file pattern for haplotype files,
##' when '.legend.gz' is added as a suffix, and `chrom` is substituted, should resolve to a
##' haplotype legend file.  Similarly when `.hap.gz` is added as a suffix should resolve to a
##' haplotype file
##' @return a dataframe or list of dataframes (see `bind_rows`)
##' @export
read_hap_blocks <- function(hap_dir,
                            ldmr,
                            hap_file_pattern = "EUR.chr{chrom}_1kg_geno") {

    all_files <- fs::dir_ls(hap_dir)
    all_f <- fs::path_file(all_files)
    ldmr <- sort(as_ldmap_region(ldmr))
    chrom <- unique(chromosomes(ldmr))
    ldmr_l <- split(ldmr, chromosomes(ldmr))
    legend_cols <- readr::cols(
                              id = readr::col_character(),
                              position = readr::col_integer(),
                              a0 = readr::col_character(),
                              a1 = readr::col_character()
                          )
    callback_factory <- function(chr, ldm) {
        function(x, pos) {
            dplyr::transmute(x,
                             id = id,
                             snp = new_ldmap_snp(chr, position, a0, a1),
                             snp_id = seq(from = pos,
                                          length.out = dplyr::n())) %>%
                dplyr::filter(snp %overlaps% ldm) %>%
                dplyr::mutate(snp = as.character(snp))
        }
    }
    map_fun <- purrr::imap_dfr
    map_fun(ldmr_l,
            function(ldmap_reg, chrom) {
                leg_file <- fs::path(hap_dir,
                                     glue::glue(hap_file_pattern),
                                     ext = "legend.gz")
                hap_file <- fs::path(hap_dir,
                                     glue::glue(hap_file_pattern),
                                     ext = "hap.gz")
                stopifnot(all(file.exists(leg_file)),
                          all(file.exists(hap_file)))

                cbf <- readr::DataFrameCallback$new(callback_factory(chrom, ldmap_reg))
                rdf <- readr::read_delim_chunked(file = leg_file,
                                                 delim = " ",
                                                 callback = cbf,
                                                 col_names = names(legend_cols$cols),
                                                 col_types = legend_cols,
                                                 skip = 1)
                dplyr::mutate(rdf, snp = as_ldmap_snp(snp),
                              hap = read_hap(hap_file,
                                             subset = snp_id)) %>%
                    dplyr::select(-snp_id)
            })
}

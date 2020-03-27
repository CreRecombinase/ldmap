chromosome_levels <- function(as_ucsc = TRUE) {
  paste0(
    ifelse(as_ucsc, "chr", ""),
    c(
      as.character(1:22),
      c(
        "X",
        "Y"
      )
    )
  )
}



##'
##' column specifier for chromosme
##'
##' @param prefix_chr boolean indicating whether or not there is a `chr` prefix
##' @param ... currently unused
##' @return `col` specification (see `?readr::col`)
##' @export
col_chromosome <- function(prefix_chr = TRUE, ...) {
  if (prefix_chr) {
    return(readr::col_factor(levels = chromosome_levels(TRUE)))
  }
  return(readr::col_factor(levels = chromosome_levels(FALSE)))
}




##' Default columns for a `bed` (genomic regions) file
##'
##' @export
bed_region_cols <- function(chrom = col_chromosome(),
                            start = readr::col_integer(),
                            end = readr::col_integer(),
                            ...) {
  argl <- rlang::list2(...)
  arg_l <- rlang::list2(
    chrom = chrom,
    start = start,
    end = end,
    !!!argl
  )
  return(rlang::exec(readr::cols, !!!arg_l))
}




##'
##' Read a bed (genomic interval, not binary ped) file
##'
##'
##' @param file path to plink bim file, `file` can be any type you can pass to
##' `readr::read_delim`
##' @param compact a boolean indicating whether to compact chrom,pos,ref and alt
##' as `ldmap_snp` (default to TRUE)
##' @param cols column specification
##' @param read_fun function to use to read data(must accept `col_names` and `col_types` arguments)
##' @param ...
##'
##' @return a `tibble` with the contents of the `bed` file
##' @export
read_bed <- function(file, compact = TRUE, cols = bed_region_cols(), read_fun = readr::read_tsv, ...) {
  bcn <- names(cols$cols)
  ret_df <- read_fun(file,
    col_names = bcn,
    col_types = cols, skip = 1L
  )
  if (compact) {
    return(compact_ldmap_region(ret_df, chrom = bcn[1], start = bcn[2], end = bcn[3]))
  }
  return(ret_df)
}


#' Write a dataframe with an ldmap_region column as a bed file
#'
#' @param df dataframe with at least one ldmap_region column
#' @param path output file to write to
#' @param region_col name of ldmap_region column
#' @param write_fun function to use for writing bed file (defautls to `readr`'s `write_tsv`)
#' @param ... arguments passed to `write_fun`
#'
#' @return
#' @export
#'
#' @examples
#'
#' temp_df <- tibble::tibble(ld = ldetect_EUR)
#' pth <- fs::file_temp(ext = "bed")
#' write_bed(temp_df, pth)
write_bed <- function(df, path, region_col = region_cols(df), write_fun = readr::write_tsv, ...) {
  write_fun(explode_ldmap_region(df, ldmap_region = region_col), path = path, ...)
}


#' Read HDF5 file with snp info
#'
#' @param h5file
#' @param ldmr ldmap_region scalar or vector
#' @param datapath datapath contain snp dataframe
#' @param snp_id boolean indicating if the row number should be returned (
#' @param ... other arguments passed to read_df_h5
#'
#' Requires `EigenH5`
#'
#' @return a dataframe with only SNPs from inside `ldmr`
#' @export
read_snp_region_h5 <- function(h5file, ldmr, datapath = "snp", snp_id=FALSE, ...) {
  if (!requireNamespace("EigenH5", quietly = TRUE)) {
      stop("Package \"EigenH5\" needed for this function to work",
           ", please install it ",
           "`(remotes::install_github(\"CreRecombinase/EigenH5\")`)", call. = FALSE)
  }
  argl <- rlang::list2(...)
  if ("ldmap_region_offset" %in% EigenH5::ls_h5(h5file)) {
    reg_df <- EigenH5::read_df_h5(h5file, "ldmap_region_offset") %>%
      dplyr::rename(mldmr = ldmr) %>%
      dplyr::filter(mldmr %overlaps% ldmr)
    if (nrow(reg_df) == 0) {
        snp_df <- rlang::exec(EigenH5::read_df_h5,
                              filename = h5file,
                              datapath = datapath,
                              subset = 1L,
                              !!!argl) %>%
        dplyr::slice(integer())
      return(snp_df)
    }
    ss <- purrr::flatten_int(map2(reg_df$offset, reg_df$datasize), function(offset, datasize) {
        seq(from = offset, length.out = datasize) + 1
    })
    snp_df <- rlang::exec(EigenH5::read_df_h5,
                          filename = h5file, datapath = datapath, subset = ss, !!!argl)
    if (snp_id) {
        if (! "snp_id" %in% colnames(snp_df)) {
            snp_df$snp_id <- ss
        }
    }
    sc <- snp_df[[snp_cols(snp_df)]]
    return(dplyr::filter(snp_df, sc %overlaps% ldmr))
  }
  if ("chrom_offset" %in% EigenH5::ls_h5(h5file)) {
    chroms <- unique(chromosomes(ldmr))
    co_df <- EigenH5::read_df_h5(h5file, "chrom_offset", subset = chroms)
    if (nrow(co_df) == 1) {
      ss <- seq(from = co_df$offset, length.out = co_df$datasize) + 1
      snp_df <- rlang::exec(EigenH5::read_df_h5, filename = h5file, datapath = datapath, subset = ss, !!!argl)
      sc <- snp_df[[snp_cols(snp_df)]]
      return(dplyr::filter(snp_df, sc %overlaps% ldmr))
    }
    else {
      return(purrr::pmap_dfr(co_df, function(offset, datasize, ...) {
        ss <- seq(from = offset, length.out = datasize)
        snp_df <- rlang::exec(EigenH5::read_df_h5, filename = h5file, datapath = datapath, subset = ss, !!!argl)
        sc <- snp_df[[snp_cols(snp_df)]]
        return(dplyr::filter(snp_df, sc %overlaps% ldmr))
      }))
    }
  }
  sinfo <- EigenH5::info_h5(h5file, datapaths = EigenH5::ls_h5(h5file, datapath, full_names = TRUE), attr = TRUE)
  snp_col <- dplyr::filter(sinfo, lengths(attributes) == 1) %>%
    dplyr::filter(purrr::map_lgl(attributes, function(x) {
      xc <- x[["class"]]
      if (is.null(xc)) {
        return(FALSE)
      }
      return("ldmap_snp" %in% xc)
    })) %>%
    dplyr::pull(name)
  stopifnot(length(snp_col) == 1)
  srv <- EigenH5::read_vector_h5(h5file, snp_col)
  sro <- srv %overlaps% sort(ldmr)
  snp_lgl <- which(sro)
  argl$subset <- snp_lgl
  retdf <- rlang::exec(EigenH5::read_df_h5, filename = h5file, datapath = datapath, !!!argl)
  return(retdf)
}

##' @title read a snp dataframe stored in HDF5 with a snp index
##'
##'
##' @param h5file hdf5 file
##' @param snplist list of rsids to read (alternatively a list of ldmap_snps
##' @param datapath hdf5 group that contains the dataframe
##' @param snp_id boolean for whether rownames are returned
##' @param rsidpath dataset name of the snp column to match
##' @param allow_missing whether missing data is allowed
##' @param ...
##' @return
##' @export
read_snp_rsid_h5 <- function(h5file, snplist, datapath = "snp", snp_id=FALSE, rsidpath="rsid", allow_missing = FALSE, ...) {
    if (!requireNamespace("EigenH5", quietly = TRUE)) {
        stop("Package \"EigenH5\" needed for this function to work",
             ", please install it ",
             "`(remotes::install_github(\"CreRecombinase/EigenH5\")`)", call. = FALSE)
    }
    argl <- rlang::list2(...)
    if(inherits(snplist,"ldmap_snp")){
        retdf <- read_snp_region_h5(h5file,
                                    ldmr = as_ldmap_region(snplist),
                                    datapath = datapath,
                                    snp_id = snp_id,
                                    ...)
        if(!allow_missing){
            sc <- retdf[[snp_cols(retdf)]]
            stopifnot(all(snplist %in% sc))
        }
        return(retdf)
    }
    snplist <- rsid2int(snplist)
    if(! rsidpath %in% EigenH5::ls_h5(h5file,datapath))
        stop("rsidpath: ",
             rsidpath,
             " not in datapath: ",
             datapath,
             "(",paste0(EigenH5::ls_h5(h5file,datapath),collapse=","),
             ")")
    rsidvec <- rsid2int(EigenH5::read_vector_h5(h5file, as.character(fs::path(datapath,rsidpath))))
    match_v <- rsidvec %in% snplist
    match_i <- which(match_v)
    if(!allow_missing){
        stopifnot(all(match_v))
    }
    retdf <- rlang::exec(EigenH5::read_df_h5,
                         filename = h5file,
                         datapath = datapath,
                         subset = match_i,
                         !!!argl)

    if(snp_id && !("snp_id" %in% colnames(retdf)))
        retdf$snp_id <- match_i
    return(retdf)
}

##' Read the contents of a fasta file as a vector of ldmap_alleles
##'
##' @param fastafile path to a fasta file (can be gzipped)
##' @param subset positions to subset in the file.  Default is the whole sequence
##' @param line_length the number of nucleotide characters per line (default is 50)
##' @return a vector ldmap_allele with the sequence contents of the file
##' @export
read_fasta <- function(fastafile, subset = NULL, line_length = 50L) {
  fasta_ext <- fs::path_ext(fastafile)
  if (fasta_ext == "gz") {
    fc <- gzfile(fastafile, "rb")
  } else {
    fc <- file(fastafile, "rb")
  }
  hdr <- readLines(fc, n = 1)
  if (!is.null(subset)) {
    file_subset <- subset + as.integer(subset / (line_length + 1))
    max_sub <- max(file_subset)
    mdata <- readBin(fc, what = raw(), n = max_sub)[file_subset]
    na <- new_ldmap_allele(mdata)
    close(fc)
    return(na)
  } else {
    if (fasta_ext == "gz") {
      nfc <- file(fastafile, "rb")
      tfp <- seek(nfc, where = -4, origin = "end")
      nub <- readBin(nfc, what = integer(), size = 4) - (nchar(hdr) + 1)
      close(nfc)
    } else {
      nub <- as.integer(fs::file_size(fastafile)) - (nchar(hdr) + 1)
    }
    file_subset <- seq_len(nub)[-seq(
      from = (line_length + 1),
      to = nub,
      by = (line_length + 1)
    )]
    retv <- new_ldmap_allele(readBin(fc, what = raw(), n = nub - 1)[file_subset])
    close(fc)
    return(retv)
  }
}

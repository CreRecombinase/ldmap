##' add jitter to cumulative genetic map values so that the results are stricrly sorted
##' @title Jitter genetic map
##' @param map cumulative genetic map values.  must by sorted (but need not be strictly sorted)
##' @return numeric vector of length `p` with genetic map values
##' @author Nicholas Knoblauch
jitter_map <- function(map,min_diff = 10 * .Machine$double.eps) {
    stopifnot(!is.unsorted(map))
    if (!is.unsorted(map, strictly = TRUE))
        return(map)
    md <- cumsum(rep(min_diff,length(map)))
    rmd <- map + md
    if (is.unsorted(rmd, strictly = TRUE)) {
        warning("map is still not strictly sorted, trying `jitter_map` with a larger `min_diff`")
        return(jitter_map(map, min_diff = min_diff * 10))
    }

    return(rmd)
}

##' @title Read 1000 genomes genetic map file
##'
##'
##' @param map_files path to one or more genetic map files
##' @param chromosomes chromosome corresponding to each genetic map file (default is seq_along(map_files))
##' @return dataframe with column `snp` with ldmap_snp (no alleles) and `map` with cumulative genetic map values
##' @export
read_1kg_maps <- function(map_files,chromosomes=seq_along(map_files)){
    stopifnot(length(map_files) == length(chromosomes))
    rc <- readr::cols(
                     position = readr::col_double(),
                     rate = readr::col_skip(),
                     map = readr::col_double(),
                     filtered = readr::col_double()
                 )
    purrr::map2_dfr(map_files, chromosomes, function(tmap,chr) {
        tmap_df <- vroom::vroom(tmap,
                                col_names = c("position", "rate", "map", "filtered"),
                                col_types = rc,
                                skip = 1) %>%
            dplyr::filter(filtered != 1) %>%
            dplyr::transmute(snp = new_ldmap_snp(chrom = chr,
                                                 pos = as.integer(position)),
                             map = map)
    })

}

##' @title Download genetic map data from the 1000 genomes
##'
##'
##' @param pop the 1000 genomes population
##' @param dest_dir directory to download the data to
##' @param kg_url url with 1000 genomes data
##' @return vector of paths files (see `read_1kg_maps`)
##' @export
dl_1kg_map <- function(pop="CEU",
                       dest_dir= tempdir(),
                       kg_url = "ftp://ftp-trace.ncbi.nlm.nih.gov/1000genomes/ftp/" ){
    if(! pop %in% names(population_Ne))
        stop("population must be one of",paste0(names(population_Ne),collapse=","))
    if(!fs::dir_exists(fs::path(dest_dir,pop))){
        if(!requireNamespace("rvest", quietly = TRUE)){
            stop("Package \"rvest\" is needed for dl_1kg_map to work.`)",call.=FALSE)
        }
        if(!requireNamespace("xml2", quietly = TRUE) ){
            stop("Package \"xml2\" is needed for dl_1kg_map to work.`)",call.=FALSE)
        }
        ftp_url <- paste0(kg_url,"technical/working/20130507_omni_recombination_rates/")
        data_page <- xml2::read_html(ftp_url) %>%
            rvest::html_nodes("body") %>%
            rvest::html_text() %>%
            readr::read_delim(delim = " ", col_names = c("priv", "oth", "ftp", "owner", "size", "Month", "day", "year", "file"))
        tar_files <- dplyr::filter(data_page,!stringr::str_starts(file,"README")) %>% dplyr::pull(file)
        names(tar_files) <- str_replace(tar_files,"^([A-Z]+)_.+","\\1")
        dl_urls <- paste0(ftp_url,tar_files)

        pops <- data_page

        tfile <- tar_files[pop]
        dest_file <- fs::path(dest_dir, tfile)
        download.file(paste0(ftp_url, tfile), destfile = dest_file)
        td <- tempdir()
        untar(dest_file, exdir = dest_dir)
    }
    fs::dir_ls(fs::path(dest_dir, pop))
    fs::path(dest_dir, pop, glue::glue("{pop}-{1:22}-final.txt.gz"))

}

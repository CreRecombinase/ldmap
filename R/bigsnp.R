


##' estimate reference panel LD
##'
##' @param gwas_df dataframe with a ldmap_snp column (and optionally an effect-size column (beta-hat)
##' @param reference_file `rds` file with `bigsnp` data
##' @param LDshrink boolean whether to use LDshrink
##' @param drop_missing whether to drop missing data or throw an error
##' @return correlation matrix
##' @author Nicholas Knoblauch
reference_panel_ld <- function(gwas_df, reference_file, LDshrink = TRUE, drop_missing = TRUE) {

    stopifnot(is.character(reference_file))
    bs <- bigsnpr::snp_attach(reference_file)
    bsx <- bs$genotypes
    bsmap <- tibble::as_tibble(bs$map) %>%
        compact_snp_struct(chrom =  "chromosome",
                           pos =  "physical.pos",
                           ref = "allele2",
                           alt = "allele1")

    ret_df <- match_ref_panel(gwas_df, bsmap$snp_struct) %>%
        dplyr::select(-match_type)
    if (!drop_missing) {
        stopifnot(all(!is.na(ret_df$index)))
    }else {
        ret_df <- dplyr::filter(ret_df, !is.na(index))
    }
    nbsx <- bsx[, ret_df$index]
    if (LDshrink) {
        ret_df$map <- bsmap$genetic.dist[ret_df$index]
        if (is.unsorted(ret_df$map, strictly = TRUE)) {
            ret_df$map <- jitter_map(ret_df$map)
        }
        stopifnot(!is.unsorted(ret_df$map, strictly = TRUE))

        R <- ldshrink::ldshrink(nbsx,ret_df$map)
    }else {
        R <- cor(nbsx)
    }
    rownames(R) <- as.character(ret_df$match)
    colnames(R) <- as.character(ret_df$match)
    return(R)
}


#' Align to bigsnp reference
#'
#' @param gwas_df gwas summary stats
#' @param reference_file bigsnp reference file
#' @param remove_missing whether to remove missing data
#' @param read_map_fun function for reading metadata from reference_file
#' @return gwas summary statistics (ready for fine-mapping etc)
#' @export
#'
align_reference <- function(gwas_df, reference_file, remove_missing = TRUE, read_map_fn =  identity){
    bsmap <- read_map_fn(reference_file)
    ## bs <- bigsnpr::snp_attach(reference_file)
    ## bsmap <- tibble::as_tibble(bs$map) %>%
    ##   compact_snp_struct(chrom =  "chromosome",
    ##                      pos =  "physical.pos",
    ##                      ref = "allele2",
    ##                      alt = "allele1",
    ##                      remove = FALSE)
    sc_df <- snp_cols(gwas_df)
    ret_df <- match_ref_panel(gwas_df,bsmap$snp_struct) %>%
        dplyr::select({{sc_df}} := match)
    if(remove_missing){
        ret_df <- dplyr::filter(ret_df,!is.na(index))
    }else{
        stopifnot(all(!is.na(ret_df$index)))
    }
    return(ret_df)

}

##' subset bigsnp by an ldmap_region
##'
##' @param ldmr an ldmap_region
##' @param reference_files reference files
##' @param output_file output file
##' @param init_fn function that chooses which of the refernce files to
##' read from reference_files must take args ldmr and reference_files
##' @param filter_map_fn function that takes one refernce file and an ldmap range and returns a subset map
##' @param filter_geno_fun funcion that takes one reference file, a dataframe with snp metadata and returns a matrix of genotypes
##' @param write_fun function for writing the data
##' @return returns results of write_fun
##' @author Nicholas Knoblauch
##' @export
subset_rds <- function(ldmr, reference_files, output_file,init_fn , filter_map_fn, filter_geno_fn, write_fn) {
    stopifnot(length(ldmr) == 1)
    reference_files <- init_fn(reference_files = reference_files, ldmr = ldmr)
    map <- filter_map_fn(reference_file = reference_files, ldmr = ldmr)
    bsx <- filter_geno_fn(reference_file = reference_files, ldmr = ldmr, map = map)
    write_fn(map, bsx, output_file)
}


##' estimate LD from a reference panel
##'
##'
##' @param reference_file path to reference haplotype data
##' @param LDshrink boolean for whether to use ldshrink
##' @param read_map_fun function for reading snp metadata
##' @param read_dosage_fun function for reading dosage data
##' @param reference_files a reference file
##' @return LD matrix
##' @author Nicholas Knoblauch
##' @export
panel_ld <- function(reference_file, LDshrink = TRUE, read_map_fn, read_dosage_fn) {

    bsmap <- read_map_fn(reference_file)
    nbsx <- read_dosage_fn(reference_file, bsmap)
    ## nbsx <- bs$genotypes[,]
    if (LDshrink) {
        if (is.unsorted(bsmap$map, strictly = TRUE)) {
            bsmap$map <- jitter_map(bsmap$map)
        }
        stopifnot(!is.unsorted(bsmap$map, strictly = TRUE))
        R <- ldshrink::ldshrink(nbsx, bsmap$map)
    }else {
        R <- cor(nbsx)
    }
    rownames(R) <- as.character(bsmap$snp_struct)
    colnames(R) <- as.character(bsmap$snp_struct)
    return(R)
}

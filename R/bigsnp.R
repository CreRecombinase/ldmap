##' estimate reference panel LD
##'
##' @param gwas_df dataframe with a ldmap_snp column (and optionally an effect-size column (beta-hat)
##' @param reference_file `rds` file with `bigsnp` data
##' @param LDshrink boolean whether to use LDshrink
##' @param drop_missing whether to drop missing data or throw an error
##' @return correlation matrix
##' @author Nicholas Knoblauch
reference_panel_ld <- function(gwas_df, reference_file,LDshrink = TRUE,drop_missing = TRUE) {
    ## library(ldmap)
    ## data(ldetect_EUR)
    ## reference_file <- example_bigsnp()
    ## drop_missing <- TRUE
    ## LDshrink <- TRUE
    ## gwas_df <- readRDS(system.file("test_gwas_df.RDS",package = "ldmap")) %>% dplyr::mutate(ldmr=snp_in_range(snp_struct,ldetect_EUR)) %>%
    ##     dplyr::filter(ldmr==unique(ldmr)[3])

    stopifnot(is.character(reference_file))
    bs <- bigsnpr::snp_attach(reference_file)
    bsx <- bs$genotypes
    bsmap <- tibble::as_tibble(bs$map) %>%
        compact_snp_struct(chrom =  "chromosome",
                           pos =  "physical.pos",
                           ref = "allele2",
                           alt = "allele1")

    ret_df <- match_ref_panel(gwas_df, bsmap$snp_struct) %>% dplyr::select(-match_type)
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

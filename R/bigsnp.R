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



subset.bigSNP2 <- function(x,
                          ind.row,
                          ind.col,
                          backingfile = bigsnpr:::getNewFile(x, "sub"),
                          ...) {

  G <- x$genotypes
  # Support for negative indices
  ind.row <- bigstatsr::rows_along(G)[ind.row]
  ind.col <- bigstatsr::cols_along(G)[ind.col]



  # Create new FBM and fill it
  G2 <- bigstatsr::FBM.code256(
    nrow = length(ind.row),
    ncol = length(ind.col),
    code = G$code256,
    init = NULL,
    backingfile = backingfile,
    create_bk = TRUE
  )
  bigsnpr:::replaceSNP(G2, G, rowInd = ind.row, colInd = ind.col)

  # http://stackoverflow.com/q/19565621/6103040
  newfam <- x$fam[ind.row, , drop = FALSE]
  rownames(newfam) <- bigstatsr::rows_along(newfam)
  newmap <- x$map[ind.col, , drop = FALSE]
  rownames(newmap) <- bigstatsr::rows_along(newmap)

  # Create the bigSNP object
  snp.list <- structure(list(genotypes = G2,
                             fam = newfam,
                             map = newmap),
                        class = "bigSNP")

  # save it and return the path of the saved object
  rds <- bigstatsr::sub_bk(G2$backingfile, ".rds")
  saveRDS(snp.list, rds)
  rds
}



#' Align to bigsnp reference
#'
#' @param gwas_df gwas summary stats
#' @param reference_file bigsnp reference file
#'
#' @return gwas summary statistics (ready for fine-mapping etc)
#' @export
#'
align_reference <- function(gwas_df,reference_file,remove_missing = TRUE){
  bs <- bigsnpr::snp_attach(reference_file)
  bsmap <- tibble::as_tibble(bs$map) %>%
    compact_snp_struct(chrom =  "chromosome",
                       pos =  "physical.pos",
                       ref = "allele2",
                       alt = "allele1",
                       remove = FALSE)
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
##' @param reference_files bigsnpr rds files
##' @return a new rds file
##' @author Nicholas Knoblauch
##' @export
subset_rds <- function(ldmr, reference_files, pattern="sub") {
    stopifnot(length(ldmr) == 1)

    bs <- bigsnpr::snp_attach(reference_files[[chromosomes(ldmr)]])
    stopifnot(all(bs$chromosome == chromosomes(ldmr)))
    bsmap <- tibble::as_tibble(bs$map) %>%
        dplyr::mutate(index = 1:dplyr::n()) %>%
        dplyr::filter(
                   dplyr::between(
                              physical.pos,
                              starts(ldmr) - 1L,
                              (ends(ldmr) - 1L))) %>%
        compact_snp_struct(chrom =  "chromosome",
                           pos =  "physical.pos",
                           ref = "allele2",
                           alt = "allele1",
                           remove = FALSE)
    bsx <- bs$genotypes

    rdsfile <- subset.bigSNP2(bs,
                              ind.row = seq_len(nrow(bsx)),
                              ind.col = bsmap$index,
                              backingfile = bigsnpr:::getNewFile(bs, pattern))
    return(rdsfile)
}


##' estimate LD from a reference panel
##'
##'
##' @param reference_files a vector of 22 bigsnpr files
##' @param LDshrink boolean for whether to use ldshrink
##' @return LD matrix
##' @author Nicholas Knoblauch
##' @export
panel_ld <- function(reference_file, LDshrink = TRUE) {

    bs <- bigsnpr::snp_attach(reference_file)
    bsmap <- tibble::as_tibble(bs$map) %>%
        compact_snp_struct(chrom =  "chromosome",
                           pos =  "physical.pos",
                           ref = "allele2",
                           alt = "allele1",
                           remove = FALSE)
    nbsx <- bs$genotypes[,]
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



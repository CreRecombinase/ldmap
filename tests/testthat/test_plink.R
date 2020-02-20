context("test plink reader/writer")


test_that("can read plink genotype data", {
    library(ldmap)
    plinkf <- fs::path_package("plink.bed", package = "ldmap")
    plink_df <- read_plink_bed(plinkf)
    gt <- plink_df$genotypes
    tgt <- gt[[1]]
    tgt[1:3]
})

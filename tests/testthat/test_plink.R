context("test plink reader/writer")


test_that("can read plink genotype data", {
    plinkf <- fs::path_package("plink.bed", package = "ldmap")
    famf <- fs::path_package("plink.fam", package = "ldmap")
    bimf <- fs::path_package("plink.bim", package = "ldmap")
    bim_df <- read_plink_bim(bimf)
    fam_df <- readr::read_delim(famf,
                                delim = " ",
                                col_names = names(fam_cols()$cols),
                                col_types = fam_cols())
    expect_equal(nrow(fam_df), 6)
    plink_df <- read_plink_bed(plinkf)
    fc <- file(plinkf, open = "rb", raw = TRUE)
    dv <- readBin(plinkf, what = raw(), 9)
    nv <- readBin(fc, what = raw(), n = 3)
    expect_equal(nv, as.raw(c(0x6c, 0x1b, 0x01)))
    nnv <- readBin(fc, what = raw(), n = 2)
    data_dv <- dv[-c(1:3)]
    expect_equal(nnv, data_dv[1:2])
    nnnv <- readBin(fc, what = raw(), n = 2)
    gt <- plink_df$genotypes
    tgt <- gt[[1]]
    expect_equal(data_dv[1:2], vctrs::vec_data(tgt))
    expect_equal(data_dv[3:4], vctrs::vec_data(gt[[2]]))
    expect_equal(data_dv[5:6], vctrs::vec_data(gt[[3]]))
})


test_that("can compute allele frequencies", {
    plinkf <- fs::path_package("plink.bed", package = "ldmap")
    plink_df <- read_plink_bed(plinkf)
    tgt <- plink_df$genotypes[[1]]
    (af_1 <- calc_AF(plink_df$genotypes[[1]]))
    expect_true(is.na(af_1))
    af_1<- calc_AF(plink_df$genotypes[[1]],TRUE)
    tgtd <- gt2double(tgt)
    calc_AF(tgtd,TRUE)
    expect_equal(af_1,    calc_AF(tgtd,TRUE))
    afv <- calc_AF(plink_df$genotypes,TRUE)
    c_afv <- purrr::map_dbl(plink_df$genotypes, ~calc_AF(gt2double(.x), TRUE))
    expect_equal(c_afv,afv)
})

test_that("can compute LD", {

    bigplink <- "/run/media/nwknoblauch/c67d20d8-a8bf-4ab9-91af-978fa311014a/1kg/1000G_EUR_Phase3_plink/1000G.EUR.QC.21.bed"
    plink_df <- read_plink_bed(bigplink)
    ld_df <- read.table("../../inst/test_ld.tsv.gz", header = TRUE, stringsAsFactors = FALSE) %>%
        tibble::as_tibble()

    fsnp <- dplyr::filter(ld_df,SNP_A %in% unique(SNP_A)[1:2])
    arsid <- unique(c(ld_df$SNP_A,ld_df$SNP_B))
    sub_pdf <- dplyr::filter(plink_df,rsid %in% arsid)

    geno2l <- purrr::map(sub_pdf$genotypes,gt2double)
    combo_df <- dplyr::select(sub_pdf,SNP_A=rsid,genotypes_A=genotypes) %>%
        dplyr::mutate(nav=NA_integer_ ) %>%
        dplyr::inner_join(dplyr::select(sub_pdf,SNP_B=rsid,genotypes_B=genotypes) %>%
                          dplyr::mutate(nav=NA_integer_),by=c("nav")) %>%
        dplyr::select(-nav) %>%
        dplyr::mutate(geno_a=purrr::map(genotypes_A,gt2double),
                      geno_b=purrr::map(genotypes_B,gt2double)) %>%
        dplyr::mutate(R=purrr::map2_dbl(geno_a,geno_b,cor)) %>%
        dplyr::left_join(ld_df,by=c("SNP_A","SNP_B")) %>%
        dplyr::select(SNP_A,SNP_B,R.x,R.y) %>%
        dplyr::filter(SNP_A!=SNP_B,!is.na(R.y))
    expect_equivalent(combo_df$R.x,combo_df$R.y,tolerance=1e-6)
})

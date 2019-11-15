context("sets")


testthat::test_that("we can make and print ldmap_ranges",{
  
  new_ldmap_range(chrom = 3L,start=4L,end = 5L)
  input_f <- fs::path_package("test_data/fourier_ls-all.bed.gz",package = "ldmap")
  ld_df <- readr::read_tsv(input_f,col_types=readr::cols(
    chr = readr::col_character(),
    start = readr::col_integer(),
    stop = readr::col_integer()
  ))
  
  big_range <- new_ldmap_range(chrom=factor(ld_df$chr),start=ld_df$start,end=ld_df$stop)
  cdf <- ldmap_range_2_data_frame(big_range)
  expect_equal(cdf$start,ld_df$start)
  expect_equal(cdf$end,ld_df$end)
  expect_equal(cdf$chrom,as.integer(factor(ld_df$chr)))
})


testthat::test_that("we can estimate LD scores",{
  skip_on_cran()
  fer <- purrr::safely(fs::path_package)("1kg_eur.tar.bz2",package="ldmap")
  if(is.null(fer$result))
    skip("skip ldscore estimation")
  untar()
  
  
  
})

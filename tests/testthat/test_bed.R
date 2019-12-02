context("sets")


testthat::test_that("we can make and print ldmap_ranges",{
  
  
  input_f <- fs::path_package("test_data/fourier_ls-all.bed.gz",package = "ldmap")
  ld_df <- read_bed(input_f,compact = FALSE)
  
  big_range <- new_ldmap_range(ld_df$chrom,start=ld_df$start,end=ld_df$end)
  expect_equal(chromosomes(big_range),as.integer(ld_df$chrom))
  cdf <- ldmap_range_2_data_frame(big_range)
  expect_equal(cdf$start,ld_df$start)
  expect_equal(cdf$end,ld_df$end)
  expect_equal(forcats::fct_drop(ld_df$chrom),cdf$chrom)
})


# 
# testthat::test_that("we can estimate LD scores",{
#   skip_on_cran()
# 
#   fldf <- fs::path_package("22.l2.ldscore.gz",package="ldmap")
#   l2c <- readr::cols(
#   chrom = readr::col_integer(),
#   rsid = readr::col_character(),
#   pos = readr::col_double(),
#   l2 = readr::col_double()
# )
#   l2df <- readr::read_tsv(fldf,col_names = names(l2c$cols),col_types = l2c,skip = 1L)
#     fer <- purrr::safely(fs::path_package)("1kg_eur.tar.bz2",package="ldmap")
#    if(is.null(fer$result))
#     skip("skip ldscore estimation")
#   td <- tempdir()
#   ret <- untar(fer$result,exdir=td)
#   files <- fs::dir_ls(fs::path(td,"1kg_eur"))
#   bedf <- fs::path(td,"1kg_eur","22")
#   plinkd <- plink2R::read_plink(root=bedf)
#   plinkf <- magrittr::set_colnames(plinkd$bim,c("chrom","rsid","map","pos","ref","alt")) %>% 
#     dplyr::as_tibble() %>% 
#     compact_snp_struct() %>% 
#     dplyr::mutate(snp_id=1:dplyr::n()) %>% 
#     dplyr::arrange(rank.ldmap_snp(snp_struct)) %>% dplyr::semi_join(l2df)
#   
#   window_r <- window_ldmap_range(ldmap_snp = plinkf$snp_struct,cm=plinkf$map,window = 2.0)
#   split_dfl <- match_ranges_snps(plinkf,window_r)
#   
#   retdf <-   purrr::map_dfr(split_dfl,function(df,N){
#     if(nrow(df)==1){
#       R <- matrix(1.0,1,1)
#     }else{
#       R <- cor(plinkd$bed[,df$snp_id])
#      }
#     l2 <- estimate_LDscores(R=R,N=N)
#     
#     return(dplyr::mutate(df,l2=
#   },N=nrow(plinkd$bed))
#   
#   
#   
#   
#   
# })

testthat::test_that("we can subset ldmap_ranges with ldmap_ranges",{
  data(hg19_sizes)
  data("ldetect_EUR")
  map_ldetect <- range_in_range(ldetect_EUR,hg19_sizes)
  checkr <- chromosomes(ldetect_EUR)
  expect_equal(map_ldetect,checkr)
})

testthat::test_that("we can subset bigsnp files with ldmap_ranges",{
  
  bs <- example_bigsnp()
  data(ldetect_EUR)
  ld22 <- ldetect_EUR[chromosomes(ldetect_EUR)==22]
  bsl <- c(purrr::map(1:21,~NULL),bs)
  sld <- sample(ld22,size=1)
  srds <- ldmap:::subset_rds(ldmr = sld,reference_files = bsl,pattern = "example_")
  sa <- bigsnpr::snp_attach(srds)
  expect_gte(min(sa$map$physical.pos),starts(sld))
  expect_lt(max(sa$map$physical.pos),ends(sld))
  tR <- ldmap::panel_ld(srds,FALSE)
})


testthat::test_that("we can find windows",{

    p <- 500
    tsnps <- sort(new_ldmap_snp(chrom = rep(1L,p),
                           pos =  sample(1L:2^28L, p, replace = F),
                           NA2N = TRUE))
    map <- cumsum(runif(p))
    window_r <- window_ldmap_range(tsnps,map)
    snp_l <- match_ranges_snps(window_r,tsnps)
    sir  <- snp_in_range(tsnps,window_r)
    df <- tibble::tibble(snp = tsnps,map = map,region = window_r,member = sir)
    cdiff <- dplyr::group_by(df, member) %>% dplyr::summarise(diff = max(map) - min(map)) %>% dplyr::pull(diff)

    })

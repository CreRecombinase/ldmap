context("sets")


testthat::test_that("distance works for ranges",{

  tsnps_a <- rsnp_region(ldetect_EUR[10],n=100,sort=TRUE)


  tsnps_b <- set_positions(tsnps_a,positions(tsnps_a)+200L)
  cdist_f <- distance(tsnps_b,tsnps_a)
  cdist_r <- distance(tsnps_a,tsnps_b)
  
  tldmr <- new_ldmap_region(chromosomes(tsnps_a),positions(tsnps_a),positions(tsnps_b)+1)
  ret <- nearest_snp_region(tsnps_a,tldmr)
  
  expect_equal(rep(200L,100),cdist_f)
  expect_equal(rep(-200L,100),cdist_r)
  
  tregions <- rregion(n=100,chroms=4L)
  tsnps <- rsnp_region(hg19_sizes[4],n=100)

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
#   window_r <- window_ldmap_region(ldmap_snp = plinkf$snp_struct,cm=plinkf$map,window = 2.0)
#   split_dfl <- match_regions_snps(plinkf,window_r)
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

testthat::test_that("we can subset ldmap_regions with ldmap_regions",{
  data(hg19_sizes)
  data("ldetect_EUR")
  adf <- ldmap_region_2_data_frame(ldetect_EUR)
  cr <- region_in_region(ldetect_EUR[250],hg19_sizes)
  map_ldetect <- region_in_region(ldetect_EUR,hg19_sizes)

  adf <- dplyr::mutate(adf,ldmr=map_ldetect,ldid=1:dplyr::n())
  checkr <- chromosomes(ldetect_EUR)
  expect_equal(map_ldetect,checkr)
})



testthat::test_that("we can find windows",{

    p <- 500
    tsnps <- sort(new_ldmap_snp(chrom = rep(1L,p),
                           pos =  sample(as.integer(2^28), p, replace = F),
                          ))
    map <- cumsum(runif(p))
    window_r <- window_ldmap_region(tsnps,map,window = 1)
    expect_true(all(starts(window_r)<=positions(tsnps)))
    expect_true(all(ends(window_r)>=positions(tsnps)))
    sir  <- snp_in_region(tsnps,window_r)
    df <- tibble::tibble(snp = tsnps,map = map,region = window_r,member = sir)
    cdiff_df <- dplyr::group_by(df, member) %>% 
      dplyr::summarise(diff = max(map) - min(map)) 
    cdiff <- cdiff_df    %>%
      dplyr::pull(diff)
    expect_true(all(cdiff<0.5))

    })

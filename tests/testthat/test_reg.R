context("assigning to ld region")

 input_f <- system.file("test_data/fourier_ls-all.bed.gz",package = "ldshrink")
  ld_df <- readr::read_tsv(input_f,col_types=readr::cols(
    chr = readr::col_character(),
    start = readr::col_integer(),
    stop = readr::col_integer()
  ))
  ld_df <- dplyr::group_by(ld_df,chr)
  ld_df <- 
    dplyr::mutate(ld_df,start = dplyr::if_else(start==min(start),
                                               0L,start),
                  stop = dplyr::if_else(stop==max(stop),
                                        .Machine$integer.max,stop))
    ld_df <- dplyr::ungroup(ld_df)
    ld_df <- dplyr::mutate(ld_df,region_id=1:dplyr::n())
  
  n_ld_df <- nrow(ld_df)
  ld_df <-   dplyr::mutate(ld_df,snp_ct=sample(1:100,n_ld_df,replace=T))
  
  snp_df <- purrr::pmap_dfr(ld_df,function(chr,start,stop,region_id,snp_ct,...){
    n_snps <- snp_ct
    dplyr::mutate(tibble::tibble(pos=sort(sample(start:stop,n_snps,replace=T))),
                  chr=chr,region_id=region_id)
  })

test_that("Check for assigning SNPs to blocks",{
  
 

  n_region_id <- assign_region(break_chr = ld_df$chr,
                               break_start = ld_df$start,
                               break_stop = ld_df$stop,
                               break_id = ld_df$region_id,
                               snp_chr = snp_df$chr,
                               snp_pos = snp_df$pos,assign_all=T)
                               
                               
  expect_equal(paste0(snp_df$region_id,".0"),n_region_id)
})


# 
# test_that("Check for finding SNPs works with all snps",{
#   
# 
#    snp_df <- purrr::pmap_dfr(ld_df,function(chr,start,stop,region_id,snp_ct,...){
#     n_snps <- snp_ct
#     dplyr::mutate(tibble::tibble(pos=sort(sample(start:stop,n_snps,replace=T))),
#                   chr=chr,region_id=region_id)
#   })
#   
#   snp_sample <- sort(sample(1:nrow(snp_df),100,replace=F))
#   query_slice <- dplyr::slice(snp_df,snp_sample)
# #   ret <- ldshrink::find_alleles(query_chrom = query_slice$chr,query_pos = query_slice$pos,ref_chrom = snp_df$chr,ref_pos = snp_df$pos,query_chunk = query_slice$region_id,ref_chunk = snp_df$region_id)
# # expect_equal(ret,snp_sample)
# })



test_that("Check for finding SNPs works with all snps with count limit",{

  max_size <- 10L
  snp_df <- dplyr::ungroup(dplyr::mutate(dplyr::group_by(snp_df,region_id),
                                           ict=1:dplyr::n(),
                                           n_reg_id=paste0(region_id,".",as.integer(ict/max_size))))
  
  
  n_region_id <- assign_region(break_chr = ld_df$chr,
                               break_start = ld_df$start,
                               break_stop = ld_df$stop,
                               break_id = ld_df$region_id,
                               snp_chr = snp_df$chr,max_size = 10L,
                               snp_pos = snp_df$pos,assign_all=T)
  #snp_df <- dplyr::mutate(snp_df,n_region_id=n_region_id,snp_id=1:dplyr::n())
                               
  expect_equal(snp_df$n_reg_id,n_region_id)
})


test_that("Check for finding SNPs works with all snps with count limit",{

  
  max_size <- 10L
  min_size <- 5L
  snp_df <- dplyr::ungroup(dplyr::mutate(dplyr::group_by(snp_df,region_id),
                                           ict=1:dplyr::n(),
                                          sreg_id= as.integer(ict/max_size),
                                           n_reg_id=paste0(region_id,".",sreg_id)))
  snp_df <- dplyr::inner_join(snp_df,dplyr::summarise(dplyr::group_by(snp_df,n_reg_id),nct=dplyr::n()))
  snp_df <- dplyr::mutate(snp_df,n_reg_id=dplyr::if_else( nct < 5L, paste0(region_id,".",pmax(sreg_id-1,0)),n_reg_id))
  n_region_id <- assign_region(break_chr = ld_df$chr,
                               break_start = ld_df$start,
                               break_stop = ld_df$stop,
                               break_id = ld_df$region_id,
                               snp_chr = snp_df$chr,max_size = 10L,min_size = 5,
                               snp_pos = snp_df$pos,assign_all=T)
  expect_equal(snp_df$n_reg_id,n_region_id)
})


test_that("Check for finding SNPs works with missing snps",{
  
  input_f <- system.file("test_data/fourier_ls-all.bed.gz",package = "ldshrink")
  ld_df <- readr::read_tsv(input_f,col_types=readr::cols(
    chr = readr::col_character(),
    start = readr::col_integer(),
    stop = readr::col_integer()
  ))
  ld_df <- dplyr::group_by(ld_df,chr)
  ld_df <- 
    dplyr::mutate(ld_df,start = dplyr::if_else(start==min(start),
                                               0L,start),
                  stop = dplyr::if_else(stop==max(stop),
                                        .Machine$integer.max,stop))
  ld_df <- dplyr::ungroup(ld_df)
  ld_df <- dplyr::mutate(ld_df,region_id=1:dplyr::n(),chr=as.integer(gsub("chr","",chr)))
  
  
  snp_df <- purrr::pmap_dfr(ld_df,function(chr,start,stop,region_id,...){
    n_snps <- sample(1:100,1)
    dplyr::mutate(tibble::tibble(pos=sort(sample(start:stop,n_snps,replace=T))),
                  chr=chr,region_id=region_id)
  })
  
  snp_df <- dplyr::mutate(snp_df,snp_id=1:dplyr::n())
  snp_sample <- sort(sample(1:nrow(snp_df),100,replace=F))
  osnp_sample <- sample(snp_sample,50,replace=F)
  query_slice <- dplyr::slice(snp_df,snp_sample)
  snp_df <- snp_df <- snp_df[-osnp_sample,]
  ret <- ldshrink::find_alleles(query_chrom = query_slice$chr,query_pos = query_slice$pos,ref_chrom = snp_df$chr,ref_pos = snp_df$pos,query_chunk = query_slice$region_id,ref_chunk = snp_df$region_id)
  nret <- dplyr::inner_join(query_slice,snp_df,by = c("pos", "chr", "region_id", "snp_id"))
  expect_equal(nrow(nret),50L)
  
})

context("assigning to ld region")


library(ldmap)

data("ldetect_EUR")
ld_df <- explode_ldmap_range(tibble::tibble(ldmr=ldetect_EUR),remove = FALSE) %>% dplyr::mutate(region_id=1:dplyr::n())
ld_df <-   dplyr::mutate(ld_df,snp_ct=sample(1:10,dplyr::n(),replace=T))

snp_df <- purrr::pmap_dfr(ld_df,function(chrom,start,end,region_id,snp_ct,...){
  n_snps <- snp_ct
  retsnps <- as.integer(sort(sample(start:(end-1L),n_snps,replace=F)))
  tibble::tibble(pos=retsnps,
                 chrom=rep(chrom,n_snps),region_id=rep(region_id,n_snps))
})

snp_df <- dplyr::mutate(snp_df,chrom=as.integer(gsub("chr","",chrom)),pos=pos) %>% compact_snp_struct(ref=NA,alt=NA,remove = FALSE)

  # test_that("New check for assignment",{
  #   
  #   snp_df <- rsnp_region(ldmap_region = head(ldetect_EUR),n = 5) %>% 
  #     tibble::enframe(name = NULL,value="snp") %>% dplyr::mutate(region_id=rep(1L:6L,each=5))
  #   
  #   id_assignment <- snp_df$region_id
  #   
  #   snp_df <- dplyr::mutate(snp_df,result=snp_in_range(ldmap_snp = snp,ldmap_range = head(ldetect_EUR))) %>% 
  #     dplyr::arrange(snp) %>% 
  #     dplyr::mutate(ldmap=cumsum(runif(dplyr::n())))
  #   tmap <- window_ldmap_range(ldmap_snp = snp_df$snp,cm = snp_df$ldmap)
  # 
  #   snp_l <- match_ranges_snps(snp_df,ld_df$ldmr)
  #   expect_equal(as_ldmap_range(names(snp_l)),ld_df$ldmr)
  #   
  #   
  #   
  #   testthat::expect_equal(unname(purrr::map_int(snp_l,nrow)),ld_df$snp_ct)
  #   testthat::expect_equal(snp_df$result,snp_df$region_id)
  #   })
  # 
  
  test_that("we can get something looking like torus-formatted input from snps and annotations",{
    
     
    k <- 10
    kv <- paste0("feature_",1:k)
    ld_df <- dplyr::mutate(ld_df,list_l=sample(kv,dplyr::n(),replace=TRUE))
    ldl <- split(ld_df$ldmr,ld_df$list_l)
    ldl <- purrr::map(ldl,sort)
    snp_df <- dplyr::inner_join(snp_df,dplyr::select(ld_df,region_id,list_l)) %>% dplyr::distinct(snp_struct,.keep_all=TRUE)
    snp_dfr <- dplyr::mutate(snp_df,value=1L) %>% dplyr::select(snp_struct,list_l,value) %>% 
      tidyr::spread(key=list_l,value=value,fill=0L) %>%
      dplyr::rename(ldmap_snp=snp_struct)
    rl <- snp_in_ranges(snp_df$snp_struct,ldl)
    expect_equal(rl,snp_dfr)
    anno_k <- 15
    rn <- sample(50:100,anno_k,replace=FALSE)
    anno_l <- purrr::map(rn,~rregion(n = .x,sort = TRUE))
    expect_equal(lengths(anno_l),rn)
    
    
    
    
  })
  
  test_that("we can concatenate vectors of ranges",{
    data("ldetect_EUR")
    ldra <- ldetect_EUR[1]
    ldrb <- ldetect_EUR[2]
    expect_equal(c(ldra,ldrb),ldetect_EUR[1:2])
    
  })
  
#   
# test_that("Check for assigning SNPs to blocks",{
#   
#   n_region_id <- assign_region(break_chr = ld_df$chrom,
#                                break_start = ld_df$start,
#                                break_stop = ld_df$stop,
#                                break_id = ld_df$region_id,
#                                snp_chr = snp_df$chr,
#                                snp_pos = snp_df$pos,assign_all=T)
#                                
#                                
#   expect_equal(paste0(snp_df$region_id,".0"),n_region_id)
# })


test_that("can merge regions",{
  
  split_ldf_1 <- ld_df$ldmr[1:nrow(ld_df)%%2==0]
  split_ldf_2 <- ld_df$ldmr[1:nrow(ld_df)%%2!=0]
  result <- merge_ldmap_ranges(split_ldf_1,split_ldf_2)  
  expect_equal(result,ld_df$ldmr)
  
})



test_that("we can correctly assign regions to regions when they're length 1",{

  for(i in seq_along(hg19_sizes)){  
    rir <- range_in_range(ldetect_EUR,hg19_sizes[i])
    expect_equal(chromosomes(ldetect_EUR)==chromosomes(hg19_sizes[i]),!is.na(rir))
  }
})
 


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


# 
# test_that("Check for finding SNPs works with all snps with count limit",{
# 
#   max_size <- 10L
#   snp_df <- dplyr::ungroup(dplyr::mutate(dplyr::group_by(snp_df,region_id),
#                                            ict=1:dplyr::n(),
#                                            n_reg_id=paste0(region_id,".",as.integer(ict/max_size))))
#   
#   
#   n_region_id <- assign_region(break_chr = ld_df$chr,
#                                break_start = ld_df$start,
#                                break_stop = ld_df$stop,
#                                break_id = ld_df$region_id,
#                                snp_chr = snp_df$chr,max_size = 10L,
#                                snp_pos = snp_df$pos,assign_all=T)
#   #snp_df <- dplyr::mutate(snp_df,n_region_id=n_region_id,snp_id=1:dplyr::n())
#                                
#   expect_equal(snp_df$n_reg_id,n_region_id)
# })


# test_that("Check for finding SNPs works with all snps with count limit",{
# 
#   
#   max_size <- 10L
#   min_size <- 5L
#   snp_df <- dplyr::ungroup(dplyr::mutate(dplyr::group_by(snp_df,region_id),
#                                            ict=1:dplyr::n(),
#                                           sreg_id= as.integer(ict/max_size),
#                                            n_reg_id=paste0(region_id,".",sreg_id)))
#   snp_df <- dplyr::inner_join(snp_df,dplyr::summarise(dplyr::group_by(snp_df,n_reg_id),nct=dplyr::n()))
#   snp_df <- dplyr::mutate(snp_df,n_reg_id=dplyr::if_else( nct < 5L, paste0(region_id,".",pmax(sreg_id-1,0)),n_reg_id))
#   n_region_id <- assign_region(break_chr = ld_df$chr,
#                                break_start = ld_df$start,
#                                break_stop = ld_df$stop,
#                                break_id = ld_df$region_id,
#                                snp_chr = snp_df$chr,max_size = 10L,min_size = 5,
#                                snp_pos = snp_df$pos,assign_all=T)
#   expect_equal(snp_df$n_reg_id,n_region_id)
# })

# 
# test_that("Check for finding SNPs works with missing snps",{
#   input_f <- fs::path_package("test_data/fourier_ls-all.bed.gz",package = "ldmap")
#   ld_df <- readr::read_tsv(input_f,col_types=readr::cols(
#     chr = readr::col_character(),
#     start = readr::col_integer(),
#     stop = readr::col_integer()
#   ))
#   ld_df <- dplyr::group_by(ld_df,chr)
#   ld_df <- 
#     dplyr::mutate(ld_df,start = dplyr::if_else(start==min(start),
#                                                0L,start),
#                   stop = dplyr::if_else(stop==max(stop),
#                                         .Machine$integer.max,stop))
#   ld_df <- dplyr::ungroup(ld_df)
#   ld_df <- dplyr::mutate(ld_df,region_id=1:dplyr::n(),chr=as.integer(gsub("chr","",chr)))
#   
#   
#   snp_df <- purrr::pmap_dfr(ld_df,function(chr,start,stop,region_id,...){
#     n_snps <- sample(1:100,1)
#     dplyr::mutate(tibble::tibble(pos=sort(sample(start:stop,n_snps,replace=T))),
#                   chr=chr,region_id=region_id)
#   })
#   
#   snp_df <- dplyr::mutate(snp_df,snp_id=1:dplyr::n())
#   snp_sample <- sort(sample(1:nrow(snp_df),100,replace=F))
#   osnp_sample <- sample(snp_sample,50,replace=F)
#   query_slice <- dplyr::slice(snp_df,snp_sample)
#   snp_df <- snp_df <- snp_df[-osnp_sample,]
#   ret <- ldshrink::find_alleles(query_chrom = query_slice$chr,query_pos = query_slice$pos,ref_chrom = snp_df$chr,ref_pos = snp_df$pos,query_chunk = query_slice$region_id,ref_chunk = snp_df$region_id)
#   nret <- dplyr::inner_join(query_slice,snp_df,by = c("pos", "chr", "region_id", "snp_id"))
#   expect_equal(nrow(nret),50L)
#   
# })

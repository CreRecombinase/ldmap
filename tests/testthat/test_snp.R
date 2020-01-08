context ("bitpacking")

testthat::test_that("combining alleles works",{
  expect_equal(
    ldmap::make_ambig("T","A"),
    ldmap::make_ambig("A","T"))
  expect_equal(
    ldmap::make_ambig("T","T"),"T")
  
  expect_equal(
    ldmap::make_ambig(
      ldmap::make_ambig(
        ldmap::make_ambig("A","C"),"T"),
      "G"),
    "N")
  
  expect_equal(extract_alt("A",make_ambig("A","T")),"T")
})

testthat::test_that("crazy bit packing works",{
  p <- 1e5
  nucs <- c("A","C","T","G","N")
  names(nucs) <- nucs
  letter_to_int <- new_ldmap_allele(1L:15L)
  chrom <- sample(1:23,p,replace=TRUE)
  pos <- as.integer(sample(2^28,p,replace=T))
  ref <- sample(letter_to_int,p,replace=T)
  alt <- sample(letter_to_int,p,replace=T)
  ret <- new_ldmap_snp(chrom,pos,ref,alt)        
  so_ret <- sort(ret)
  cdf <- ldmap::ldmap_snp_2_dataframe(ret)
  testthat::expect_equal(cdf$chrom,chrom)
  testthat::expect_equal(cdf$pos,pos)
  testthat::expect_equal(cdf$ref,ref)
  testthat::expect_equal(cdf$alt,alt)
})



testthat::test_that("sorting works",{
  p <- 1e5
  nucs <- c("A","C","T","G","N")

  names(nucs) <- nucs
  letter_to_int <- new_ldmap_allele(nucs)
  chrom <- sample(1:23,p,replace=TRUE)
  pos <- sample(as.integer(max(ends(hg19_sizes))),p,replace=T)
  ref <- sample(letter_to_int,p,replace=T)
  alt <- sample(letter_to_int,p,replace=T)
  ret <- new_ldmap_snp(chrom,pos,ref,alt)      
  hd <- as.double(ret)
  hd <- vctrs::vec_cast(ret,double())
  shd <- sort(hd)
  tret <- as_ldmap_snp(hd)
  ldshd <- as_ldmap_snp(shd)
  so_ret <- sort(ret)
  vso_ret <- vctrs::vec_sort(ret)
  expect_true(!is.unsorted(chromosomes(so_ret)))
  
  expect_equal(tret,ret)
  
  expect_equal(hd,unclass(ret))
  
  ntib <- tibble::tibble(id=ret,chrom,pos,ref,alt) %>% 
    dplyr::mutate(order_ret=order.ldmap_snp(id))
  sort_ntib <- dplyr::arrange(ntib,chrom,pos,ref,alt)
  expect_equal(so_ret,sort_ntib$id)
})





testthat::test_that("matching works like bigsnpstatsr",{
  

  # gwasf <- "/home/nwknoblauch/Dropbox/scratch/ptb_scratch/ptb_gwas.h5"
  # codf <- EigenH5::read_df_h5(gwasf,"chrom_offset")
  # 
  # tdf <- EigenH5::read_df_h5(gwasf,"snp",offset=14379541L,datasize=170459L) %>% dplyr::mutate(
  #   snp_struct=new_ldmap_snp(chrom,
  #                            pos,
  #                            ascii_ref = ref,
  #                            ascii_alt = alt
  #                            )) %>% dplyr::select(snp_struct,beta,pval) %>% dplyr::arrange(rank.ldmap_snp(snp_struct))
  # 
  # 
  tdf <- readRDS(fs::path_package("test_gwas_df.RDS",package = "ldmap")) %>% tibble::as_tibble()
  
  geno_df <-   readRDS(fs::path_package("test_reference.RDS",package = "ldmap")) %>% 
    tibble::as_tibble() %>% explode_snp_struct(alleles_to_character = TRUE,remove=FALSE)
    
  
  sumstats <- explode_snp_struct(tdf,alleles_to_character = TRUE) %>% 
    dplyr::select(chr=chrom,
                            pos,a0=ascii_ref,a1=ascii_alt,beta,p=pval) %>% dplyr::mutate(snp_id=1:dplyr::n())
  
  rs <-  readRDS(fs::path_package("snp_match.RDS",package = "ldmap")) %>% tibble::as_tibble()
  
  
  
  info_snp_id <- geno_df$snp_struct
  
  

  jm <- join_snp(query = tdf$snp_struct,
                 reference = info_snp_id,
                 rsid = as.integer(gsub(pattern = "rs",replacement = "",x = geno_df$id))) %>% 
    dplyr::mutate(snp_id=1:dplyr::n())
  mjm <- dplyr::filter(jm,!is.na(index)) %>% dplyr::select(snp_struct=match,rsid)
  testthat::expect_equal(mjm,dplyr::select(rs,snp_struct,rsid))
})




testthat::test_that("we can index with ourselves",{
  p <- 1e3
  nucs <- c("A","C","T","G","N")
  names(nucs) <- nucs
  letter_to_int <- new_ldmap_allele(nucs)
  chrom <- sample(1:23,p,replace=TRUE)
  pos <- sample(2^43,p,replace=T)
  ref <- sample(letter_to_int,p,replace=T)
  alt <- sample(letter_to_int,p,replace=T)
  ret <- new_ldmap_snp(chrom,pos,ref,alt)        
  sret <- sort(ret[-1])
  sret_id <- 1:length(sret)
  query <- sort(c(sample(ret[-1],size=100),ret[1]))
  
  perf_match <- join_snp(query,sret)
  
  perf_match_rsid <- join_snp(query,sret,rsid = sret_id)
  spmatch_rsid <- perf_match_rsid[!is.na(perf_match_rsid$index),]
  expect_equal(sret[spmatch_rsid$index],spmatch_rsid$match)
  
  expect_equal(nrow(perf_match),101)
  expect_equal(query[perf_match$match_type=="snp_match"],
               perf_match$match[perf_match$match_type=="snp_match"],check.attributes = FALSE)
  
  
  
  
})

testthat::test_that("we get the best match possible",{
  p <- 1e3
  nucs <- c("A","C","T","G","N")
  names(nucs) <- nucs
  letter_to_int <- new_ldmap_allele(nucs)
  chrom <- sample(1:23,1)
  pos <- sample(2^43,1)
  ref <- new_ldmap_allele("G")
  alt <- new_ldmap_allele("A")
  ref_q <- new_ldmap_snp(chrom,pos,ref,alt)
  pmatch <- new_ldmap_snp(chrom,pos-1,alt,ref) #doesn't match because it goes before
  gmatch <- new_ldmap_snp(chrom,pos,alt,ref) #reverse match
  bmatch <- new_ldmap_snp(chrom,pos,new_ldmap_allele("T"),new_ldmap_allele("C")) #matches on position and not on alleles
  amatch <- new_ldmap_snp(chrom,pos+1,alt,ref) #doesn't match because it goes after
  
  matches <- c(pmatch,gmatch,bmatch,amatch)
  matches2 <- c(pmatch,bmatch,gmatch,amatch)
  result <- join_snp(ref_q,matches)
  rev_result <- join_snp(ref_q,matches2[1:2])
  expect_equal(result$match,rev_result$match)
  
})

test_that("we can migrate old SNP storage to new SNP storage",{
  fldf <- fs::path_package("22.l2.ldscore.gz",package="ldmap")
  l2c <- readr::cols(
    chrom = readr::col_integer(),
    rsid = readr::col_character(),
    pos = readr::col_double(),
    l2 = readr::col_double()
  )
  l2df <- readr::read_tsv(fldf,col_names = names(l2c$cols),col_types = l2c,skip = 1L)
  tsnps <- new_ldmap_snp(chrom=l2df$chrom,pos=l2df$pos)
  osnps <- readRDS(fs::path_package(package = "ldmap","old_snps.RDS"))
  csnps <- ldmap:::migrate_ldmap_snp(osnps)
  
  expect_equal(chromosomes(csnps),chromosomes(tsnps))
  expect_equal(positions(csnps),positions(tsnps))
})




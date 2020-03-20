## code to prepare `ld_df` dataset goes here
#create a dataframe 
                                        #


input_f <- fs::path_package("test_data/fourier_ls-all.bed.gz",package = "ldmap")
ld_df <- read_bed(input_f,compact = FALSE)
ld_df <- dplyr::group_by(ld_df,chrom)
ld_df <- dplyr::mutate(ld_df,
                       start = dplyr::if_else(start==min(start),
                                              1L,start),
                       stop = dplyr::if_else(end==max(end),
                                             as.integer(2^29-1),end)) %>%
    dplyr::ungroup() %>% dplyr::arrange(chrom,start,end)
  
ldetect_EUR <- ldmap::new_ldmap_region(ld_df$chrom, start=ld_df$start, end=ld_df$end)
usethis::use_data(ldetect_EUR, overwrite = TRUE)

chrom_df <- readr::read_tsv("http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes",col_names=c("chrom","end")) %>% 
  dplyr::filter(chrom %in% paste0("chr",c(1:22,"X"))) %>% 
  dplyr::mutate(chrom=factor(chrom,levels=paste0("chr",c(1:22,c("X","Y"))))) %>% dplyr::arrange(chrom)
hg19_sizes <- new_ldmap_region(chrom=chrom_df$chrom,start=rep(1L,nrow(chrom_df)),end=chrom_df$end+1L)
usethis::use_data(hg19_sizes,overwrite = TRUE)

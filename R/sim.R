#' Generate random genomic regions
#'
#' @param n number of ranges to generate
#' @param chroms chromosomes from which to sample
#' @param sort whether to sort the output
#'
#' @return a vector of length `n` of (optionally sorted ldmap_regions)
#' @export
#'
#' @examples
#' #chromosomes are assumed to be hg19
#' head(rregion(100,chroms=1L))
rregion <- function(n,chroms=1L:23L,sort=TRUE){
  stopifnot(all(chroms %in% 1L:23L))
  hg19_df <- ldmap_region_2_data_frame(hg19_sizes)
  if(length(chroms)==1){
    chrom <- rep(chroms,n)
  }else{
    chrom <- sample(chroms,n,replace=T)
  }
  
  tdf <- dplyr::mutate(tibble::tibble(chrom=chrom),
                         start=as.integer(runif(n=dplyr::n(),min=1,max=hg19_df$end[chrom])),
                         end=as.integer(runif(n=dplyr::n(),min=start+1,max = hg19_df$end[chrom])))
  if(sort)
    tdf <- dplyr::arrange(tdf,chrom,start,end)
  
  stopifnot(all(tdf$end>tdf$start))
  return(new_ldmap_region(tdf$chrom,tdf$start,tdf$end))
}
##' Draw SNPs from random positions in the genome
##'
##'
##' @param n number of SNPs to simulate
##' @param chroms chromosomes from which to draw
##' @param sort whether to sort output
##' @param ... currently unused
##' @return
##' @author Nicholas Knoblauch
rsnp <- function(n, chroms = 1L:23L, sort = TRUE,replace=FALSE, ...){
  stopifnot(all(chroms %in% 1L:23L))
  hg19_df <- ldmap_region_2_data_frame(hg19_sizes)
  hg19_df$end  <-  as.integer(hg19_df$end - 1)
  chrom <- sample(chroms, n, replace = T)
  
  tdf <- dplyr::mutate(
                    tibble::tibble(chrom = chrom),
                    pos = sample(1L:hg19_df$end, size = dplyr::n(),replace = replace))
  if (sort)
    tdf <- dplyr::arrange(tdf, chrom, pos)

  return(new_ldmap_snp(tdf$chrom, tdf$pos))
  
}


#' Draw random snps from ldmap regions
#'
#' @param ldmap_region ldmap region to draw samples from
#' @param n number of draws, recycled to match ldmap_region
#' @param sort whether to sort the result
#'
#' @return vector of ldmap_snps
#' @export
#'
#' @examples
#' #generate one snp per ldetect block
#' data(ldetect_EUR)
#' snp_per_block <- rsnp_region(ldetect_EUR,n=1)
#' ##generate three snps from one ldetect block
#' snps <- rsnp_region(ldetect_EUR[4],n=3)
#' ##generate 10 snps in total, 3 from one ldetect block
#' ##four from the second and 3 from the third
#' snps <- rsnp_region(ldetect_EUR[4:6],n=c(3,4))
rsnp_region <- function(ldmap_region, n, sort = TRUE,replace=FALSE){
  stopifnot(length(n)<=length(ldmap_region))
    n <- as.integer(n[(seq_along(ldmap_region) %% length(n))+1])
    chrom <- rep(chromosomes(ldmap_region),n)
    pos <- purrr::flatten_int(purrr::map2(ldmap_region,n,function(x,y){
      sample(as.integer(starts(x)):as.integer(ends(x)),y,replace=replace)
    }))
    rets <- new_ldmap_snp(chrom,pos)
    if(sort){
      return(sort(rets))
    }
    return(rets)
}




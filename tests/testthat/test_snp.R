context ("bitpacking")




testthat::test_that("crazy bit packing works",{
        p <- 1e5
        nucs <- c("A","C","T","G","N")
        names(nucs) <- nucs
        letter_to_int <- purrr::map_int(nucs,utf8ToInt)
        chrom <- sample(1:23,p,replace=TRUE)
        pos <- sample(2^43,p,replace=T)
        ref <- sample(letter_to_int,p,replace=T)
        alt <- sample(letter_to_int,p,replace=T)
        ret <- new_ldmap_snp(chrom,pos,ref,alt)        
        so_ret <- sort(ret)
        head(sort(ret))
        cdf <- ldmap::ldmap_snp_2_dataframe(ret)
        testthat::expect_equal(cdf$chrom,chrom)
        testthat::expect_equal(cdf$pos,pos)
        testthat::expect_equal(cdf$ascii_ref,unname(ref))
        testthat::expect_equal(cdf$ascii_alt,unname(alt))
})



testthat::test_that("sorting works",{
        p <- 1e5
        nucs <- c("A","C","T","G","N")
        names(nucs) <- nucs
        letter_to_int <- purrr::map_int(nucs,utf8ToInt)
        chrom <- sample(1:23,p,replace=TRUE)
        pos <- sample(2^43,p,replace=T)
        ref <- sample(letter_to_int,p,replace=T)
        alt <- sample(letter_to_int,p,replace=T)
        ret <- new_ldmap_snp(chrom,pos,ref,alt)        
        hd <- as.double(ret)
        hd <- vctrs::vec_cast(ret,double())
        tret <- as_ldmap_snp(hd)
        expect_equal(tret,ret)
        expect_equal(hd,unclass(ret))
        
        ntib <- tibble::tibble(id=ret,chrom,pos,ref,alt) %>% dplyr::mutate(order_ret=order.ldmap_snp(id))
        sort_ntib <- dplyr::arrange(ntib,chrom,pos,ref,alt)
        sort_ntib2 <- dplyr::slice(ntib,order(id))
        
        expect_equal(sort_ntib2,sort_ntib,ignore_row_order=FALSE)
        sort_ntib2 <- dplyr::arrange(ntib,rank.ldmap_snp(id),ref,alt)
        expect_equal(sort_ntib2,sort_ntib,ignore_row_order=FALSE)
        
        
        
        
        op <- order(chromosomes(ret),positions(ret))
        order_ret <- order.ldmap_snp(ret)
        so_ret <- sort(ret)
        
        od_ret <- ret[order.ldmap_snp(ret)]
        expect_true(!is.unsorted(chromosomes(so_ret)))
        expect_true(!is.unsorted(chromosomes(od_ret)))
        expect_true(all(split(positions(so_ret),chromosomes(so_ret)) %>% purrr::map_lgl(~!is.unsorted(.x))))
        expect_true(all(split(positions(od_ret),chromosomes(od_ret)) %>% purrr::map_lgl(~!is.unsorted(.x))))
        expect_equal(as.numeric(unclass(ret[op])),as.numeric(unclass(ret[order_ret])))
        
})




testthat::test_that("we can index with ourselves",{
        p <- 1e3
        nucs <- c("A","C","T","G","N")
        names(nucs) <- nucs
        letter_to_int <- purrr::map_int(nucs,utf8ToInt)
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
        expect_equal(sret[perf_match_rsid$rsid],perf_match_rsid$match)
        
        expect_equal(nrow(perf_match),101)
        expect_equal(sret[perf_match$index],perf_match$match)
        expect_equal(perf_match$query[perf_match$match_type=="snp_match"],
                     perf_match$match[perf_match$match_type=="snp_match"],check.attributes = FALSE)
        
        
        
        
})

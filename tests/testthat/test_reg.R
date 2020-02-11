context("assigning to ld region")


data("ldetect_EUR")
ld_df <- explode_ldmap_region(tibble::tibble(ldmr = ldetect_EUR), remove = FALSE) %>%
    dplyr::mutate(region_id = 1:dplyr::n())
ld_df <- dplyr::mutate(ld_df, snp_ct = sample(1:10, dplyr::n(), replace = T))

snp_df <- purrr::pmap_dfr(ld_df, function(chrom, start, end, region_id, snp_ct, ...) {
    n_snps <- snp_ct
    retsnps <- as.integer(sort(sample(start:(end - 1L), n_snps, replace = F)))
    tibble::tibble(
        pos = retsnps,
        chrom = rep(chrom, n_snps), region_id = rep(region_id, n_snps)
    )
})

snp_df <- dplyr::mutate(snp_df, chrom = as.integer(gsub("chr", "", chrom)), pos = pos) %>%
    compact_snp_struct(ref = NA, alt = NA, remove = FALSE)

test_that("we can get something looking like torus-formatted input from snps and annotations", {
    k <- 10
    kv <- paste0("feature_", 1:k)
    ld_df <- dplyr::mutate(ld_df, list_l = sample(kv, dplyr::n(), replace = TRUE))
    ldl <- split(ld_df$ldmr, ld_df$list_l)
    ldl <- purrr::map(ldl, sort)
    snp_df <- dplyr::inner_join(snp_df, dplyr::select(ld_df, region_id, list_l)) %>%
        dplyr::distinct(snp_struct, .keep_all = TRUE)
    snp_dfr <- dplyr::mutate(snp_df, value = 1L) %>%
        dplyr::select(snp_struct, list_l, value) %>%
        tidyr::spread(key = list_l, value = value, fill = 0L) %>%
        dplyr::rename(ldmap_snp = snp_struct)
    rl <- snp_in_regions(snp_df$snp_struct, ldl)
    expect_equal(rl, snp_dfr)
    anno_k <- 15
    rn <- sample(50:100, anno_k, replace = FALSE)
    anno_l <- purrr::map(rn, ~ rregion(n = .x, sort = TRUE))
    expect_equal(lengths(anno_l), rn)
})

test_that("we can concatenate vectors of ranges", {
    data("ldetect_EUR")
    ldra <- ldetect_EUR[1]
    ldrb <- ldetect_EUR[2]
    expect_equal(c(ldra, ldrb), ldetect_EUR[1:2])
})



test_that("can merge regions", {
    split_ldf_1 <- ld_df$ldmr[seq_len(nrow(ld_df)) %% 2 == 0]
    split_ldf_2 <- ld_df$ldmr[seq_len(nrow(ld_df)) %% 2 != 0]
    result <- merge_ldmap_regions(split_ldf_1, split_ldf_2)
    expect_equal(result, ld_df$ldmr)
})



test_that("we can bind_rows", {
    lda <- tibble::tibble(lda = ldetect_EUR[12])
    ldb <- tibble::tibble(lda = ldetect_EUR[15])

    ret_df <- purrr::map_df(ldetect_EUR[c(12, 15)], ~ tibble::tibble(lda = .x))
    ldr <- dplyr::bind_rows(lda, ldb)
    expect_equal(tibble::tibble(lda = ldetect_EUR[c(12, 15)]), ret_df)
    expect_equal(ret_df, ldr)
})


test_that("no overlap grouping is 1:n",{
    
    gr <- group_regions(ldetect_EUR)
    expect_equal(seq_along(ldetect_EUR),gr)
})

test_that("all overlap grouping is rep(1,n) (per chromosome)",{
    lde <- set_ends(ldetect_EUR,ends(ldetect_EUR)+1)
    ngr <- group_regions(lde)
    expect_equal(ngr,chromosomes(lde))
})

test_that("we can correctly assign regions to regions without overlap when length 1", {
    for (i in seq_along(hg19_sizes)) {
        rir <- region_in_region(ldetect_EUR, hg19_sizes[i])
        expect_equal(chromosomes(ldetect_EUR) == chromosomes(hg19_sizes[i]), !is.na(rir))
    }

})

test_that("we can correctly assign regions to regions without overlap", {
    for (i in seq_along(hg19_sizes)) {
        rir <- region_in_region(ldetect_EUR, hg19_sizes[i], TRUE)
        expect_equal(chromosomes(ldetect_EUR) == chromosomes(hg19_sizes[i]), !is.na(rir))
    }
})


test_that("overlap operators work", {
    all_t <- snp_df$snp_struct %overlaps% ldetect_EUR
    expect_true(all(all_t))
    for (i in seq_along(ldetect_EUR)) {
        rir <- snp_df$snp_struct %overlaps% ldetect_EUR[-i]
        expect_equal(rir, snp_df$region_id != i)
    }
})

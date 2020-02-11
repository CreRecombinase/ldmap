library(ldmap)
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

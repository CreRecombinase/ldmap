context("htslib")

test_that("we can open and close bgzf files",{

    bgzff <- "/home/nwknoblauch/tmp/age.ldsc.imputed_v3.male.tsv.bgz"
    fd <- ldmap:::open_bgzf(bgzff)
    num_b <- ldmap:::num_bgzf_blocks(fd)
    testthat::expect_equal(num_b,5864)
    rl <- purrr::flatten_chr(ldmap:::readlines_bgz(fd))
    ldmap:::format_bgzf(fd)

    rgz <- gzfile(bgzff)


    ir <- ldmap:::read_bgzf(fd)
    ldmap:::format_bgzf(fd)
    gd <- ldmap:::get_bgzf_data(fd)
    gdl <- ldmap:::readlines_chunk_bgzf(fd)
    tgd <- strsplit(gd,split="\n")[[1]]
    expect_equal(tgd,gdl)
    # ir <- ldmap:::read_bgzf(fd)
    # ngd <- ldmap:::get_bgzf_data(fd)
    # ngdl <-ldmap:::readlines_chunk_bgzf(fd)
    # expect_equal(ir, 0)


})

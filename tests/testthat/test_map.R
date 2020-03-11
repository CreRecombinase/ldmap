context("test interpolation of genetic map")


test_that("we can interpolate a single value between two values", {
    library(ldmap)
    tsnps <- new_ldmap_snp(1L, c(10L, 20L))
    maps <- c(1.0, 2.0)

    qsnp <- new_ldmap_snp(1L, c(15L))
    tres <- ldmap::new_interpolate_genetic_map(maps, tsnps, qsnp)
    testthat::expect_equal(tres, 1.5)

    qsnp <- new_ldmap_snp(1L, c(9L))
    tres <- ldmap::new_interpolate_genetic_map(maps, tsnps, qsnp)
    testthat::expect_equal(tres, .9)

    qsnp <- new_ldmap_snp(1L, c(21L))
    tres <- ldmap::new_interpolate_genetic_map(maps, tsnps, qsnp)
    testthat::expect_equal(tres, 2.1)

    qsnp <- new_ldmap_snp(1L, c(25L))
    tres <- ldmap::new_interpolate_genetic_map(maps, tsnps, qsnp)
    testthat::expect_equal(tres, 2.5)

    tsnps <- new_ldmap_snp(c(1L, 1L, 2L), c(10L, 20L, 30L))
    maps <- c(1.0, 2.0, 1.0)

    qsnp <- new_ldmap_snp(1L, c(15L))
    tres <- ldmap::new_interpolate_genetic_map(maps, tsnps, qsnp)
    testthat::expect_equal(tres, 1.5)

    qsnp <- new_ldmap_snp(1L, c(20L))
    tres <- ldmap::new_interpolate_genetic_map(maps, tsnps, qsnp)
    testthat::expect_equal(tres, 2.0)

    qsnp <- new_ldmap_snp(1L, c(21L))
    tres <- ldmap::new_interpolate_genetic_map(maps, tsnps, qsnp)
    testthat::expect_equal(tres, 2.1)
})




test_that("we can interpolate a two values between two values", {

    tsnps <- new_ldmap_snp(1L, c(10L, 20L))
    maps <- c(1.0, 2.0)

    qsnp <- new_ldmap_snp(1L, c(15L, 16L))
    tres <- ldmap::new_interpolate_genetic_map(maps, tsnps, qsnp)
    testthat::expect_equal(tres, c(1.5, 1.6))

    qsnp <- new_ldmap_snp(1L, c(9L, 21L))
    tres <- ldmap::new_interpolate_genetic_map(maps, tsnps, qsnp)
    testthat::expect_equal(tres, c(0.9, 2.1))

    qsnp <- new_ldmap_snp(1L, c(21L, 22L))
    tres <- ldmap::new_interpolate_genetic_map(maps, tsnps, qsnp)
    testthat::expect_equal(tres, c(2.1, 2.2))

    qsnp <- new_ldmap_snp(1L, c(25L, 26L))
    tres <- ldmap::new_interpolate_genetic_map(maps, tsnps, qsnp)
    testthat::expect_equal(tres, c(2.5, 2.6))

    tsnps <- new_ldmap_snp(c(1L, 1L, 2L), c(10L, 20L, 30L))
    maps <- c(1.0, 2.0, 1.0)

    qsnp <- new_ldmap_snp(1L, c(15L,16L))
    tres <- ldmap::new_interpolate_genetic_map(maps, tsnps, qsnp)
    testthat::expect_equal(tres, c(1.5,1.6))

    qsnp <- new_ldmap_snp(1L, c(20L,21L))
    tres <- ldmap::new_interpolate_genetic_map(maps, tsnps, qsnp)
    testthat::expect_equal(tres, c(2.0,2.1))

    qsnp <- new_ldmap_snp(1L, c(21L))
    tres <- ldmap::new_interpolate_genetic_map(maps, tsnps, qsnp)
    testthat::expect_equal(tres, 2.1)
})

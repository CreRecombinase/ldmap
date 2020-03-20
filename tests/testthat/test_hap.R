context("haplotype")


test_that("we can create and format ldmap_ht", {

    tx <- c(0, 1, 0, 0, 0)
    x <- new_ldmap_ht(tx)
    expect_equal(nchar(format(x)), 5)
    iss <- sample(c(0, 1), 64, replace = TRUE)
    tx <- new_ldmap_ht(iss)
    expect_equal(length(tx), 1)
    tix <- structure(vctrs::vec_data(tx), class = "integer64")
    expect_equal(bit64::as.bitstring(tix), format(tx))

})


mymap <- dl_1kg_map("CEU", dest_dir = "/run/media/nwknoblauch/Data/1kg")
mymap_df <- read_1kg_maps(mymap)
leg_df <- read_hap_legend(fs::path_package("EUR.chr19.legend.gz", package = "ldmap")) %>%
    assign_genetic_map(mymap_df)
samp_df <- read_hap_samples(fs::path_package("EUR.chr19.samples", package = "ldmap"))
hap_l <- read_hap(fs::path_package("EUR.chr19.hap", package = "ldmap"))
thap <- matrix(scan(fs::path_package("EUR.chr19.hap", package = "ldmap"),
                    what = integer(),
                    sep = " "
                    ),
               nrow = 250,
               byrow = TRUE
               )






test_that("we can read haplegend data", {
    expect_equal(nrow(leg_df), length(hap_l))
})


test_that("we can compute covariance", {
    th <- t(thap)
    tS <- cov(t(thap))
    cS <- ldmap:::cov_ht(hap_l[[1]], hap_l[[2]])
    expect_equal(tS[1, 2], cS)
    mcS <- ldmap:::cov_htm(hap_l)
    expect_equal(tS, mcS)
})


test_that("we can compute LDshrink like the toy R implementation", {

    calcLDR <- function(hmata, mapa, m = 85, Ne = 11490.672741, cutoff = 0.001, isgeno = FALSE) {
        S <- stats::cov(hmata)
        p <- length(mapa)
        td <- abs(outer(mapa, mapa, `-`))
        if (isgeno) {
            S <- 0.5 * S
        }

        rho <- 4 * Ne * (td) / 100
        rho <- -rho / (2 * m)
        tshrinkage <- exp(rho)
        tshrinkage[tshrinkage < cutoff] <- 0
                                        # diag(tshrinkage) <- 1
        S <- S * tshrinkage
        theta <- calc_theta(m)

        eye <- diag(p) * (0.5 * theta * (1 - 0.5 * theta))
        SigHat <- ((1 - theta) * (1 - theta)) * S + eye
        return(stats::cov2cor(SigHat))
    }

    tLDR <- calcLDR(t(thap), leg_df$map)
    mLD <- ldmap:::ldshrink_S(hap_l, leg_df$map)
    expect_equal(mLD, tLDR, tolerance = 1e-4)
})

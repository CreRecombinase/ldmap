##' pack an array of 1s and 0s into integers
##'
##' .. content for \details{} ..
##' @param x data consisting of  all `0`s and `1`s
##' @param compact_64 whether to return an `int64` or not
##' @return an integer array of `ceiling(length(x)/int_size)` where
##' `int_size is 32 unless `compact_64` is TRUE, in which case it's `64`
##' @author Nicholas Knoblauch
compact_haplotype <- function(x, compact_64 = FALSE) {

    stopifnot(all(x %in%  c(0L, 1L)))
    bx <- bit::as.bit(x)
    length(bx)


}

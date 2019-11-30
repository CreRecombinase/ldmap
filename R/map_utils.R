##' add jitter to cumulative genetic map values so that the results are stricrly sorted
##' @title Jitter genetic map
##' @param map cumulative genetic map values.  must by sorted (but need not be strictly sorted)
##' @return numeric vector of length `p` with genetic map values
##' @author Nicholas Knoblauch
jitter_map <- function(map,min_diff = 10 * .Machine$double.eps) {
    stopifnot(!is.unsorted(map))
    if (!is.unsorted(map, strictly = TRUE))
        return(map)
    md <- cumsum(rep(min_diff,length(map)))
    rmd <- map + md
    if (is.unsorted(rmd, strictly = TRUE)) {
        warning("map is still not strictly sorted, trying `jitter_map` with a larger `min_diff`")
        return(jitter_map(map, min_diff = min_diff * 10))
    }

    return(rmd)
}

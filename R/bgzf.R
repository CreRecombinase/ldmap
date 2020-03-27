
readlines_bgz <- function(filename){
    if(is.character(filename)){
        fd <- open_bgzf(filename)
    }else{
        stopifnot(typeof(filename) == "externalptr")
        fd <- filename
    }
    input_len <- seq_len(num_bgzf_blocks(fd))
    retl <- purrr::map(input_len,function(x){
        read_bgzf(fd)
        return(readlines_chunk_bgzf(fd))
    })
    close_bgzf(fd)
    return(retl)
}

format_spdi <- function(x,assembly="hg19"){
  if(tolower(assembly)%in%c("grch37","hg19")){
    data("GRCh37_metadata")
    df <- GRCh37_metadata
  }else{
    stopifnot(tolower(assembly)%in%c("grch38"))
    df <- GRCh38_metadata
  }
  nchroms <- df$RefSeq.Accn
  chroms <- chromosomes(x)
  retvec <- glue::glue("{nchroms[chroms]}:{positions(x)}:1:{alt_alleles(x)}")
  retvec
  
}


get_rsid <- function(spdi){
  url <- httr::modify_url("https://api.ncbi.nlm.nih.gov",path=as.character(fs::path("/variation/v0/spdi/",spdi,"rsids")))
  res <- httr::GET(url,httr::accept_json())
  httr::stop_for_status(res)
  ct <- httr::content(res)
  
}

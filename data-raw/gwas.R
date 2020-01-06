library(tidyverse)
fldf <- fs::path_package("22.l2.ldscore.gz",package="ldmap")
  fldf <- fs::path_package("22.l2.ldscore.gz",package="ldmap")
  l2c <- readr::cols(
  chrom = readr::col_integer(),
  rsid = readr::col_character(),
  pos = readr::col_double(),
  l2 = readr::col_double()
)
l2df <- readr::read_tsv(fldf,col_names = names(l2c$cols),col_types = l2c,skip = 1L)
snps <- new_ldmap_snp(chrom=l2df$chrom,pos=l2df$pos)

saveRDS(snps,"inst/old_snps.RDS")
fer <- fs::path_package("1kg_eur.tar.bz2",package="ldmap")
td <- tempdir()
ret <- untar(fer,exdir=td)
files <- fs::dir_ls(fs::path(td,"1kg_eur"))
bimf <- fs::path(td,"1kg_eur","22.bim")
bimc <- cols(
    chrom = col_integer(),
    rsid =  col_character(),
    map = col_double(),
    pos =  col_double(),
    alt = col_character(),
    ref = col_character()
)
bim_df <- read_tsv(bimf)

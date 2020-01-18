#chr_names_url <- 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_assembly_report.txt'
#
#
library(tidyverse)
chrf_38 <- "~/Downloads/GCF_000001405.39_GRCh38.p13_assembly_report.txt"
chrf_37 <- "~/Downloads/GCF_000001405.25_GRCh37.p13_assembly_report.txt"
chl_38 <- read_lines(chrf_38)
chl_37 <- read_lines(chrf_37)
header_38 <- tail(chl_38[str_detect(chl_38,"^#")],1) %>% str_remove("# ") %>% str_split(pattern="\t") %>% flatten_chr() %>% make.names()
header_37 <- tail(chl_37[str_detect(chl_37,"^#")],1) %>% str_remove("# ") %>% str_split(pattern="\t") %>% flatten_chr() %>% make.names()

GRCh37_metadata <- read_tsv(chrf_37,comment="#",col_names=header_37)
GRCh38_metadata <- read_tsv(chrf_38,comment="#",col_names=header_38)
usethis::use_data(GRCh37_metadata)
usethis::use_data(GRCh38_metadata)


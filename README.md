
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ldmap

<!-- badges: start -->

<!-- badges: end -->

The goal of ldmap is to simplify common operations on statistical
genetics and genomics datasets by providing a set of simple, performant
classes for dealing with (human) genomic annotation data and GWAS
summary statistics.

## Installation

<!-- You can install the released version of ldmap from [CRAN](https://CRAN.R-project.org) with: -->

<!-- ``` r -->

<!-- install.packages("ldmap") -->

<!-- ``` -->

You can install the development version of ldmap from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("CreRecombinase/ldmap")
```

# Some cool things you can do with ldmap

## Read and manipulate genomic intervals

```r

rawbed <- r"(track name=pairedReads description="Clone Paired Reads" useScore=1
chr22	1000	5001
chr22	2000	6000)" 
ft <- fs::file_temp()
write(rawbed, ft)
bed_data <- ldmap::read_bed(ft)
dplyr::filter(bed_data, widths(ldmap_region) == 4000)
```


## calculate LD for a region of the genome from plink files

```r

dl_url <- "https://data.broadinstitute.org/alkesgroup/LDSCORE/1kg_eur.tar.bz2" #example file with chromosome 22
dest_f <- fs::file_temp("1kg_eur", ext = ".tar.bz2")
download.file(
    dl_url,
    dest_f
) # The Price are kind enough to host the 1000G Phase3 genotypes for EUR
temp_d <- fs::path_dir(dest_f)
dest_d <- utils::untar(dest_f, exdir = fs::path_dir(dest_f))
dest_files <- fs::dir_ls(fs::path(temp_d, "1kg_eur"))
bim_file <- fs::path(temp_d, "1kg_eur","22",ext="bim")
bed_file <- fs::path_ext_set(bim_file, "bed")
bim_df <- ldmap::read_plink_bim(bim_file) # read plink bim file into a dataframe
bim_subset <- which(dplyr::between(
                               ldmap::positions(bim_df$snp_struct),
                               left = 16896762,
                               right = 18966860
))
plink_geno_l <- ldmap::read_plink_bed(bed_file, subset = bim_subset)
geno_matrix <- ldmap::gt2matrix(plink_geno_l)
R <- cor(geno_matrix)
image(R)
```



<!-- ```{r example} -->

<!-- library(ggplot2) -->

<!-- library(patchwork) -->

<!-- p1 <- ggplot(mtcars) + geom_point(aes(mpg, disp))+ylab("short label") -->

<!-- p2 <- ggplot(mtcars) + geom_boxplot(aes(gear, disp, group = gear))+ylab("really really\nreally long\nlabel ") -->

<!-- p1 +theme(axis.title.y = element_text(vjust=1,debug=TRUE))+ p2+theme(axis.title.y = element_text(vjust=0,debug=TRUE))+plot_layout(ncol=1,guides = "keep",) -->

<!-- ## basic example code -->

<!-- ``` -->

<!-- ```{r} -->

<!-- p1 + p2+plot_layout(ncol=1,guides = "keep")&theme(axis.text.y = element_text(hjust=1)) -->

<!-- ``` -->

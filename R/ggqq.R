manhattan_ggplot <- function(data_df,snp_col=snp_struct,p_col=pval){


## chrom_sizes <- mutate(data_df, chrom = chromosomes(snp_struct)) %>%
##     group_by(chrom) %>%
##     summarise(size=max(positions({{snp_col}}))) %>%
##     mutate(start_p=cumsum(lag(size+(48116032/2),default=0)))

## ggpmdf <- mdf %>%
##     left_join(top_snp_kg) %>%
##     mutate(
##         pos = positions(snp_struct)+chrom_sizes$start_p[chromosomes(snp_struct)],
##         chr_color = as.integer((chromosomes(snp_struct) %% 2) + 1)
##     ) %>% mutate(chr_color=if_else(is.na(name),chr_color,3L))
## isdf <- filter(ggpmdf, !is.na(name)) %>% mutate(name = str_remove(name, ".+_"))

## chrom_ticks <- mutate(ggpmdf, chrom = chromosomes(snp_struct)) %>%
##     group_by(chrom) %>%
##     summarise(med_pos = median(pos))

## ggplot(ggpmdf, aes(x = pos, y = -log10(pval), col = factor(chr_color))) +
##     theme_minimal(base_family = "Helvetica") +
##     geom_point() +
##     scale_color_manual(values = cchrom_cols) +
##     scale_x_continuous(name = "Chromosome",
##                        breaks = chrom_ticks$med_pos,
##                        labels = as.character(chrom_ticks$chrom)) +
##     geom_hline(yintercept = -log10(5e-8)) +
##     theme(
##         plot.margin = margin(),
##         panel.grid.major = element_blank(),
##         legend.position = "none",
##         panel.grid.minor = element_blank(),
##         axis.line = element_line(colour = "black"),
##         axis.ticks = element_line(colour = "black")
##     ) +
##     geom_text_repel(data = isdf, aes(x = pos, y = -log10(pval), label = name),col = "black", force = 9) +
##     ylab(bquote(-log[10](p)))
}

#' ---
#' title: Difflines variant calling stats
#' author: KDM
#' date: 2021-10-25
#' ---
#'
#'

library(tidyverse)
library(foreach)

refs = c('at6137', 'at6923', 'at6929', 'at7143', 'at8285',
         'at9104', 'at9336', 'at9503', 'at9578', 'at9744',
         'at9762', 'at9806', 'at9830', 'at9847', 'at9852',
         'at9879', 'at9883', 'at9900', 'TAIR10')

samples = foreach(varcall=c("freebayes", "mpileup"), .combine=rbind) %:%
foreach(ref=refs, .combine=rbind) %do% {
    path = sprintf("data/%s~bwa~%s~difflines~filtered-default_samples.tsv.xz", varcall, ref)
    read_tsv(path) %>%
        mutate(varcall=varcall, ref=ref)
}

snps = foreach(varcall=c("freebayes", "mpileup"), .combine=rbind) %:%
foreach(ref=refs, .combine=rbind) %do% {
    path = sprintf("data/%s~bwa~%s~difflines~filtered-default_snps.tsv.xz", varcall, ref)
    read_tsv(path) %>%
        mutate(varcall=varcall, ref=ref)
}

table(snps$metric)

snps %>%
    filter(metric == "miss") %>%
    ggplot(aes(x=percent, y=nsnp)) +
    geom_bar(aes(fill=varcall, colour=varcall), position="dodge",
             stat="identity") +
    facet_grid(ref~varcall, scales="free_y")

snp.sum = snps %>%
    filter(metric == "miss") %>%
    group_by(varcall, ref) %>%
    filter(percent < 10) %>%
    summarise(nsnps_90 = sum(nsnp))

snp.sum %>%
    ggplot(aes(x=ref, y=nsnps_90)) +
        geom_bar(aes(colour=varcall, fill=varcall), stat="identity", position="dodge")


samples %>%
    ggplot(aes(x=ref, y=missing_prop)) +
    geom_violin(aes(fill=varcall))


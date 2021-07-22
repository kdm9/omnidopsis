# # Kraken taxonomic profiling of omnidopsis
#
# Kraken uses kmer matching to identify each read's likely taxonomic origin.
# Here we summarise these reports across the full omnidopsis sample set as a
# measure of sample qualtiy. This is of particular importance for
# field-collected samples, which are likely to have non-trivial quantities of
# microbial contamination.

library(tidyverse)
library(data.table)
library(foreach)
library(doParallel)
registerDoParallel()

extract.kraken = function(file) {
    colnam = c("pct.reads.below", "n.reads.below", "n.reads",
               "level", "taxid",  "name")
    raw = fread(file, col.names=colnam)
    return(raw)
}

all.dat = foreach(file=Sys.glob("data/OA*.txt"),
                  .combine=bind_rows, .multicombine=T) %dopar% {
    res = extract.kraken(file)
    oaid = sub(".txt", "", basename(file))
    res$oaid = oaid
    res
}


min.pct = 1
levels.include = c("U", "P")
taxid.include = c(72274, 1385, 3701, 0)

Ath = readLines("data/Athaliana.txt")


select.dat = all.dat %>%
    filter(pct.reads.below > min.pct & level%in%levels.include,
           oaid %in% Ath) %>%
    mutate(name=fct_relevel(name, "unclassified", "Streptophyta"))

ggplot(select.dat, aes(y=pct.reads.below, x=oaid)) +
    geom_col(aes(fill=name, colour=name)) +
    scale_colour_brewer(palette="Paired") +
    scale_fill_brewer(palette="Paired") +
    coord_flip() +
    theme_bw() +
    labs(y="Percentage of reads", x=NULL, colour="Taxon", fill="Taxon")
ggsave("kraken-everything.svg", height=320, width=8, limitsize=F)


bad.samp = all.dat %>%
    filter(taxid==3699, pct.reads.below < 60, oaid %in% Ath) %>%
    pull(oaid)
length(bad.samp)

dat = all.dat %>%
    filter(
       pct.reads.below > min.pct & level%in%levels.include,
       oaid %in% bad.samp, oaid %in% Ath) %>%
    mutate(name=fct_relevel(name, "unclassified", "Streptophyta"))

ggplot(dat, aes(y=pct.reads.below, x=oaid)) +
    geom_col(aes(fill=name, colour=name)) +
    scale_colour_brewer(palette="Paired") +
    scale_fill_brewer(palette="Paired") +
    coord_flip() +
    theme_bw() +
    labs(y="Percentage of reads", x=NULL, colour="Taxon", fill="Taxon")
ggsave("kraken-bad.svg", height=26, width=8, limitsize=F)


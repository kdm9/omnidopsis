# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .R
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.11.2
# ---

# # Kraken taxonomic profiling of omnidopsis
#
# Kraken uses kmer matching to identify each read's likely taxonomic origin. Here we summarise these reports across the full omnidopsis sample set as a measure of sample qualtiy. This is of particular importance for field-collected samples, which are likely to have non-trivial quantities of microbial contamination.

# # Technical details

library(tidyverse)
library(data.table)
library(foreach)
library(doParallel)
registerDoParallel()

# ## Extract per-sample kraken results
#
# Kraken reports are delimited by one or more spaces, which the readers from
# tidyverse and base R can't do, so we use fread from data.table.

extract.kraken = function(file) {
    colnam = c("pct.reads.below", "n.reads.below", "n.reads",
               "level", "taxid",  "name")
    raw = fread(file, col.names=colnam)
    return(raw)
}


# Here we extract the per-sample kraken reports into one massive data frame
# with a sample column.

all.dat = foreach(file=Sys.glob("data/*.report"),
                  .combine=bind_rows, .multicombine=T) %dopar% {
    res = extract.kraken(file)
    oaid = sub(".report", "", basename(file))
    res$oaid = oaid
    res
}
write_tsv(all.dat, "out/kraken-results.tsv.xz")


# ## Sumarise to some relevant taxa/levels

# We show all samples, summarised at the phylum level.

min.pct = 1
levels.include = c("U", "P")
#taxid.include = c(72274, 1385, 3701, 0)
Ath = readLines("data/Athaliana.txt")


select.dat = all.dat %>%
    filter(pct.reads.below > min.pct & level%in%levels.include,
           oaid %in% Ath | grepl("^CAO", oaid)) %>%
    mutate(name=fct_relevel(name, "unclassified", "Streptophyta"))

ggplot(select.dat, aes(y=pct.reads.below, x=oaid)) +
    geom_col(aes(fill=name, colour=name)) +
    scale_colour_brewer(palette="Paired") +
    scale_fill_brewer(palette="Paired") +
    coord_flip() +
    theme_bw() +
    labs(y="Percentage of reads", x=NULL, colour="Taxon", fill="Taxon")
ggsave("out/kraken-everything.svg", height=320, width=8, limitsize=F)

# And now we make something the same as the above, but for only the lines with
# contamination.

bad.samp = all.dat %>%
    filter(taxid==3699, pct.reads.below < 60,
           oaid %in% Ath | grepl("^CAO", oaid)) %>%
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
ggsave("out/kraken-bad.svg", height=26, width=8, limitsize=F)


# # Just cao samples

meta = read_tsv("../../rawdata/metadata/omniath_all_accessions.tsv", guess_max=1e5)
csids = readLines("data/csids.txt")
cao = meta %>%
    filter(cs_number %in% csids)


dat = all.dat %>%
    filter(
       pct.reads.below > min.pct & level%in%levels.include,
       oaid %in% cao$oa_id | grepl("^CAO", oaid)) %>%
    mutate(name=fct_relevel(name, "unclassified", "Streptophyta"))

ggplot(dat, aes(y=pct.reads.below, x=oaid)) +
    geom_col(aes(fill=name, colour=name)) +
    scale_colour_brewer(palette="Paired") +
    scale_fill_brewer(palette="Paired") +
    coord_flip() +
    theme_bw() +
    labs(y="Percentage of reads", x=NULL, colour="Taxon", fill="Taxon")
ggsave("out/kraken-cao.svg", height=21, width=8, limitsize=F)

ggplot(dat, aes(y=n.reads.below, x=oaid)) +
    geom_col(aes(fill=name, colour=name)) +
    scale_colour_brewer(palette="Paired") +
    scale_fill_brewer(palette="Paired") +
    coord_flip() +
    theme_bw() +
    labs(y="N reads", x=NULL, colour="Taxon", fill="Taxon")
ggsave("out/kraken-cao-reads.svg", height=21, width=8, limitsize=F)

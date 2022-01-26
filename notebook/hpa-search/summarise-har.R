#' # Har search
#'
#' We wish to find which localities seem to have Har reads among the pooled Ath
#' data from the Roux lab's populations.


library(tidyverse)
library(ggmap)
library(ggrepel)
library(rnaturalearth)
library(rnaturalearthdata)

sf::sf_use_s2(FALSE)

#' First, read in the per-contig mapping statistics. Contigs are named like
#' species_contig so we split this into two columns.
fp = read_tsv("outputs/bwa~ath_har~FP.samtools.idxstats") %>%
    mutate(chr=ifelse(chr=="*", "unknown_unmapped", chr))  %>%
    separate(chr, into=c("spp", "chr"), sep="_")

#' Summarise to a per-sample percentage of mapped reads that are from Har vs
#' total mapped reads.

fp.samp = fp %>%
    group_by(sample, spp) %>%
    summarise(reads_mapped=sum(reads_mapped)) %>%
    pivot_wider(names_from="spp", values_from="reads_mapped") %>%
    transmute(sample, har_pct_mapped=100*HaR/(HaR+Ath))


#' ## Plot on a map
#'
#' We want to see if the abundance of Har reads follows any geographic trend.

meta = read_tsv("../../rawdata/metadata/omniath_all_accessions.tsv")

fp.samp.meta = fp.samp %>%
    left_join(meta, by=c("sample"="oa_id"))

#' First, what's the distribution of Har abundance across all sites.

hist(fp.samp$har_pct_mapped)

#' Mostly very low, with a density around a few percent and one sample with a lot of Har.

bbox = c(left=-1, top=45, right=4, bottom=42)
baselayer = get_stamenmap(bbox, zoom = 8, maptype = "terrain-background")
ggmap(baselayer)  +
    geom_point(aes(x=longitude, y=latitude, colour=har_pct_mapped), size=3,
               position="jitter", data=fp.samp.meta) +
    scale_colour_continuous(trans="log", name="Pct. Har Reads",
                            breaks=c(0.1, 0.3, 1, 3, 10, 30) ) +
    labs(x="Longitude", y="Latitude") +
    theme(legend.position="right")

#' There's no clear pattern across geography.

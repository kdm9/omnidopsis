# # Difflines all-against-all BAM stats
# 

remotes::install_github("mcjmigdal/sumsamstats")
library(foreach)
library(doParallel)
library(tidyverse)
registerDoParallel()

difflines_samp = readLines("../../rawdata/metadata/samplesets/difflines.txt")

allss = foreach(statfile=dir("data"), .combine=rbind) %dopar% {
    sn = sumsamstats::readSamtoolsStats(sprintf("data/%s", statfile), section="SN")$SN

    spl = strsplit(sub(".samtools.stats$", "", statfile), "~")

    sn %>%
        mutate(description = sub(":$", "", description) %>%
                             sub("1st", "first", .) %>%
                             sub("%", "", .) %>%
                             janitor::make_clean_names(.),
               reference=spl[[1]][2],
               value=as.numeric(value),
               sample=spl[[1]][3])
}

unique(allss$description)

more_metrics = allss %>%
    pivot_wider(names_from=description, values_from=value) %>%
    mutate(percentage_mapped_reads_mq0 = (reads_mq0/reads_mapped) * 100,
           percentage_reads_unmapped = (reads_unmapped/sequences) * 100)  %>% 
    pivot_longer(names_to="description", values_to="value", cols=!c(reference, sample))


metrics = c("average_quality",
            "error_rate",
            "percentage_mapped_reads_mq0", 
            "percentage_reads_unmapped",
            "percentage_of_properly_paired_reads",
            "insert_size_average")

plot_one = function(df, var) {
    title = gsub("_", " ", var) %>% tools::toTitleCase()
    ggplot(df, aes(sample, reference)) +
        geom_raster(aes(fill=value)) +
        theme_bw() + 
        theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
              legend.title=element_blank()) +
        labs(x=NULL, y=NULL, title=title)
}

plots = more_metrics %>%
    filter(sample %in% difflines_samp,
           description %in% metrics) %>%
    group_by(description) %>%
    nest() %>%
    mutate(plots = purrr::map2(data, description, plot_one))

cowplot::plot_grid(plotlist=plots$plots, ncol=1)

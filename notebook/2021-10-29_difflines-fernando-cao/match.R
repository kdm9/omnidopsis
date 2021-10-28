library(tidyverse)

meta = read_tsv("../../rawdata/metadata/omniath_all_accessions.tsv")

difflines_ecoid = readLines("../../rawdata/metadata/samplesets/difflines.txt") %>%
    sub("OAKG", "", .) %>%
    as.numeric()

cao_all = meta %>%
    filter(dataset %in% c("CR", "KG")) %>%
    mutate(csnum = cs_number %>%
                    sub("CS", "", .) %>%
                    as.numeric()) %>%
    filter(csnum <= 76425, csnum >= 76347) %>%
    arrange(ecotype_id) %>%
    mutate(in_difflines = ecotype_id %in% difflines_ecoid)


cao_dl = cao_all %>%
    filter(in_difflines)
    

write_csv(cao_dl, "cao_difflines.tsv")

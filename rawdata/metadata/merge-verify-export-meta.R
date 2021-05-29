# # Omni-ath metadata import, validation, and merge
#
# 

library(tidyverse)
library(sp)
#install.packages(c("tidyverse", "sp", "sf", "ggplot2", "ggmap", "validate"))


kg_sra = read_csv("source/1kg/1kg-sra-run-table.csv")
kg_acc = read_csv("source/1kg/1kg-accessions.csv")

af_sra = read_tsv("source/africans/african_sra-run-table-PRJEB19780.tsv")
af_acc = read_csv("source/africans/tabula-pnas.1616736114.sapp.csv", na='-')


# ## Argentinian lines
#

arg_sra = read_tsv("source/argentinians/SraRunTable_PRJEB9862.txt")

arg_acc = read_csv("source/argentinians/tabula-mec14107-sup-0001-supinfo.csv")

arg_acc = arg_acc %>%
    mutate(
        latitude = as.numeric(sp::char2dms(latitude, chd='d', chm='m')),
        longitude = as.numeric(sp::char2dms(longitude, chd='d', chm='m')),
        name=sprintf("Pat-%d", site)
    ) %>%
    select(name, latitude, longitude, elevation, collector)



arg_acc


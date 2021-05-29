# # Omni-ath metadata import, validation, and merge
#
# 

library(tidyverse)
#install.packages(c("tidyverse", "sp", "sf", "ggplot2", "ggmap", "validate"))


kg_sra = read_csv("source/1kg/1kg-sra-run-table.csv")
kg_acc = read_csv("source/1kg/1kg-accessions.csv")

african_sra = read_tsv("source/africans/african_sra-run-table-PRJEB19780.tsv")
african_acc = read_csv("source/africans/tabula-pnas.1616736114.sapp.csv", na='-')






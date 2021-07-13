# # Post-analysis of metadata

all_acc = read_tsv("omniath_all_accessions.tsv", guess_max=10000)

# # Mapping of accessions

globe_bbox = c(left=-170, bottom=-58, right=179.99, top=75)
baselayer = get_stamenmap(globe_bbox, zoom = 3, maptype = "toner-background")

ggmap(baselayer) +
    geom_point(aes(x=longitude, y=latitude, colour=dataset), data=all_acc) +
    labs(x="Longitude", y="Latitude") +
    theme(legend.position="right")
ggsave("all_accessions.svg", width=16, height=14, units="in", dpi=600)

ath_acc = all_acc %>%
    filter(species=="Arabidopsis thaliana")
ggmap(baselayer) +
    geom_point(aes(x=longitude, y=latitude, colour=dataset), data=ath_acc) +
    labs(x="Longitude", y="Latitude") +
    theme(legend.position="right")
ggsave("ath_accessions.svg", width=16, height=14, units="in", dpi=600)

# ## Tabulation

all_acc %>%
    mutate(species=sub(" subsp. .*", "", species)) %>%
    group_by(species) %>%
    summarise(n=n()) %>%
    arrange(-n) %>%
    knitr::kable()

# # SRA size and coverage estimation


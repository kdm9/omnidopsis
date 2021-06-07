# # Omni-ath metadata import, validation, and merge
#
# 
# ## OA ids
#
# Through this whole excercise, we will encounter many weird and wonderful ways
# of referring to accessions: ecotype ids, names of a million different
# schemes, and the SRA IDs of each sequencing run (which aren't 1:1 with
# accessions anyway). To smooth over this mess, I'm asigning an "oa" ID to each
# accession. These IDs consist of three parts: "oa", then a two letter dataset
# code, then a zero-padded sequential ID within that dataset. For the 1kg data,
# the numeric ID will actually just be the ecotype ID. Thus, `oakg0088` is from
# the 1kg dataset, ecotype_id 88, and `oaaf0001` is the first sample from the
# **af**rican sample set from the PNAS paper out of Angela Hancock's lab.

# ## Setup
#
# Import packages, declare helper functions

library(tidyverse)
library(readxl)
library(janitor)
library(ggmap)
library(validate)
library(sp)
library(SRAdb)
library(parzer)
#BiocManager::install("SRAdb")
#install.packages(c("tidyverse", "sp", "sf", "ggplot2", "ggmap", "validate",
#                   "readxl", "janitor"))


dist_to_land_km = function(eastings, northings, crs=4326) {
    world = sf::st_as_sf(rworldmap::getMap(resolution="high"))
    everything = data.frame(eastings, northings, i=seq_along(eastings))
    points = everything %>%
        filter(!is.na(eastings), !is.na(northings)) %>%
        sf::st_as_sf(coords=1:2, crs=crs)
    points$nearest = sf::st_nearest_feature(points, world)
    points$dists = sf::st_distance(points, world[points$nearest,], by_element=T) %>%
        units::set_units("km") %>%
        as.numeric()
    everything = dplyr::left_join(everything, points, by='i')
    everything$dists
}

get_failing_tests = function(validator_values) {
    val = values(validator_values)
    test.names = colnames(val)
    failing_tests = NULL
    val[is.na(val)] = T
    for (x in seq_len(nrow(val))) {
        failing_tests = c(failing_tests, paste(test.names[!val[x, ]], collapse=","))
    }
    failing_tests
}


multiplicity = function(x) {
    x = as.character(x)
    table(x)[x]
}

# ## Load SRA DB

sra.sql = "/home/kevin/data/work/2021-06-06_SRAmetadb.sqlite.gz"
sra = dbConnect(SQLite(), sra.sql)

# ## The OG 1001 Genomes
#
# These are the orginal 1135 genomes sequenced by a variety of labs and
# published in the 1001g consortium paper in Cell in 2016.  NB: this data
# extracted from Joffery's 1kg website, and from a SRA Run table
# next to the reads on the ebio storage unit.

kg_acc = read_csv("source/1kg/1kg-accessions.csv")
#str(kg_acc[])

kg_acc_pre = kg_acc %>%
    mutate(oa_id = sprintf("OAKG%04d", ecotype_id)) %>%
    select(oa_id, ecotype_id, sample_name=name, latitude, longitude, locality, 
           country, collector, collection_date, cs_number)
#str(kg_acc_pre)


kg_sra = read_csv("source/1kg/1kg-sra-run-table.csv")
#str(kg_sra[])

kg_sra_pre = kg_sra %>%
    select(ecotype_id=Ecotype, sra_run=Run, bioproject=BioProject,
           library_layout=LibraryLayout, instrument=Instrument,
           species=Organism)
stopifnot(all(kg_sra_pre$ecotype_id %in% kg_acc_pre$ecotype_id))
stopifnot(all(kg_acc_pre$ecotype_id %in% kg_sra_pre$ecotype_id))
#str(kg_sra_pre)


kg_meta = inner_join(kg_acc_pre, kg_sra_pre, by="ecotype_id")
str(kg_meta)

all_meta = kg_meta

# ## African lines
#
#  From Angela Hancock's lab. Durvasula & Fulgione et al. (2017)
#  doi:10.1073/pnas.1616736114 
#
# The two metadata sources aren't complete, so we will only keep the complete
# observations (inner join).

af_acc = read_csv("source/africans/tabula-pnas.1616736114.sapp.csv", na='-')
#str(af_acc[])

af_acc_pre = af_acc %>%
    mutate(collection_date = as.character(collection_year)) %>%
    select(sample_name=name, latitude=Latitude, longitude=Longitude,
           locality=Region, country=Country, collection_date) %>%
    mutate(oa_id=sprintf("OAAF%04d", 1:n()))
#str(af_acc_pre)

af_sra = read_csv("source/africans/SraRunTable.txt")
#str(af_sra[])

af_sra_pre = af_sra %>%
    select(sra_run=Run, bioproject=BioProject,
           sample_name=Alias, library_layout=LibraryLayout,
           instrument=Instrument, species=Organism)
all(af_sra_pre$sample_name %in% af_acc_pre$sample_name)
#str(af_sra_pre)

af_meta = full_join(af_acc_pre, af_sra_pre, by="sample_name")

af_missing_meta = af_meta %>%
    filter(is.na(latitude) | is.na(sra_run))
#str(af_missing_meta)

# Not all recods in the metadata were sequenced. The above is a list of the
# non-intersecting records. Keep just the complete rows (below)

af_meta_complete =  af_meta %>%
    filter(!is.na(latitude) & !is.na(sra_run))
#str(af_meta_complete)

all_meta = bind_rows(all_meta, af_meta_complete)


# ## Argentinian lines
#
# From [Luciana Kasulin's paper in MolEcol](https://dx.doi.org/10.1111/mec.14107)

arg_sra = read_tsv("source/argentinians/SraRunTable_PRJEB9862.txt")
#str(arg_sra[])

arg_sra_pre = arg_sra %>%
    select(sra_run=run_accession, bioproject=study_accession,
           sample_name=sample_alias, library_layout,
           instrument=instrument_model, species=scientific_name)
#str(arg_sra_pre)


arg_acc = read_csv("source/argentinians/tabula-mec14107-sup-0001-supinfo.csv")
#str(arg_acc[])

# The metadata need a bit of fixing. the lat and longs are in there as deg,
# min, sec, and the metadata is in there per site. we need to make it per
# accession (just add the Pat- prefix) and convert lat/long to numeric.

arg_acc_pre = arg_acc %>%
    transmute(
        sample_name=sprintf("Pat-%d", site),
        latitude = parse_lat(latitude),
        longitude = parse_lon(longitude),
        elevation,
        collector,
        country="Argentina",
        locality="Patagonia",
    )
#str(arg_acc_pre)


arg_meta = full_join(arg_sra_pre, arg_acc_pre, by="sample_name") %>%
    mutate(oa_id=sprintf("OAAR%04d", 1:n()))
str(arg_meta)


# Not all sequencing sets have metadata, so filter down to those that have
# accession metadata.

arg_meta = arg_meta %>%
    filter(!is.na(latitude))

all_meta = bind_rows(all_meta, arg_meta)


# ## Chinese lines
#
# From Yangtze river basin, see [Zou et al.
# 2017](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1378-9).

chi_sra = read_tsv("source/chinese/SraRunTable_PRJNA293798.txt")
#str(chi_sra[])
chi_sra_pre = chi_sra %>%
    select(sra_run=run_accession, bioproject=study_accession,
           sample_name=sample_alias, library_layout,
           instrument=instrument_model, species=scientific_name)
#str(chi_sra_pre)

chi_acc = read_tsv("source/chinese/supp-table-3.tsv")
#str(chi_acc[])
chi_acc_pre = chi_acc %>%
    transmute(
        sample_name=SampleID,
        latitude=Latitude,
        longitude=Longitude,
        locality=Region,
        country="China",
    )
chi_acc_pre

chi_meta = full_join(chi_acc_pre, chi_sra_pre, by="sample_name")

chi_meta %>%
    filter(is.na(latitude) | is.na(sra_run))

# This metadata is beautiful, no gaps or weird shit going on.


chi_meta = chi_meta %>%
    mutate(oa_id=sprintf("OACN%04d", 1:n()))
all_meta = bind_rows(all_meta, chi_meta)


# ## Catalonia
#
# https://www.pnas.org/content/115/52/E12443
# There are a few A. croatica, so we need to add the
# species column from SRA to the accessions table, as that doesn't have the
# species noted. 

cat_sra = read_tsv("source/catalonian/catalonianSraRunTable.txt")
#str(cat_sra[])
cat_sra_pre = cat_sra %>%
    select(sra_run=Run, bioproject=BioProject, sample_name=Sample_Name,
           library_layout=LibraryLayout, instrument=Instrument,
           species=Organism)
#str(cat_sra_pre)


cat_acc_orig = read_xlsx("source/catalonian/pnas.1816964115.sd01.xlsx") %>%
    clean_names()
#str(cat_acc_orig)

# So this metadata is a bit shitty, most of the data rows are
# actually missing metadata, and gps coords are concatentated. We
# first extract the per-deme metadata, fix it, then merge back.
#
# Also, NB: deme PO1 was in the spreadsheet incorrectly. Someone
# seems to have dragged it down, and therefore was in there as PO1,
# PO2, PO3, PO4 for the four samples. I fixed it in the spreadsheet.

cat_acc_sample = cat_acc_orig %>%
    select(sample_name, deme_of_origin)

cat_acc_deme = cat_acc_orig %>%
    select(deme_of_origin, locality=town, gps_coordinates, inland_or_coastal=location) %>%
    filter(!is.na(locality)) %>%
    unique() %>%
    mutate(gps = purrr::map(gps_coordinates, function(x) {
            g = str_split(x, ",") %>%
                unlist() %>%
                purrr::map_chr(str_trim) %>%
                as.numeric()
            data.frame(latitude=g[1], longitude=g[2])
    })) %>%
    unnest(gps) %>%
    select(-gps_coordinates)
#str(cat_acc_deme)

cat_acc = cat_acc_sample %>%
    left_join(cat_acc_deme, by="deme_of_origin")
#str(cat_acc)

cat_acc_pre = cat_acc %>%
    select(sample_name, latitude, longitude, locality)
cat_acc_pre


# What's going on here. There are extra strings in the plant IDs in the SRA
# metadata, and there are demes (pops) in the SRA that aren't in the paper's
# metadata.

# Just match names
cat_sra_pre_name = cat_sra_pre %>%
    select(sra_run, sample_name) %>%
    # nearly all samples are like deme-numericdeme-plantnum, but the JBB
    # samples are like JBBT{1,2}. Regex the JBB samples to match the rest of
    # them.
    mutate(name_fixed = ifelse(grepl("^JBBT\\d$", sample_name, perl=T),
                               sub("JBBT(\\d)", "JBB-NA-\\1", sample_name),
                               sample_name)) %>%
    # split the deme-number-plantnum into columns
    separate(name_fixed, into=c("deme", "numdeme", "plant"),
             sep="-", remove=F, extra="merge") %>%
    arrange(deme, plant) %>%
    # new name that matches the format of the table from their supps.
    mutate(match_name = sprintf("%s-%s", deme, plant))


# Merge both datasets after above "fixing". Here you can see the issue quite
# well, there are ~70 lines that match by name, but then a dozen or so that are
# either in the SRA data and not in the accession metadata, or vice versa.
# Eight of these are demes that don't have any records in the accession
# metadata (O3 and AM).

cat_merged = full_join(cat_sra_pre_name, cat_acc, by=c("match_name"="sample_name")) %>%
    arrange(is.na(deme), is.na(deme_of_origin), sample_name)

# As some of the missing data seems to be differently named plants from the
# same deme, let's merge it by deme, in essence copying the lats & longs from
# the accesssion metadata deme-wise to the SRA metadata. 

cat_acc_merge_by_deme = cat_sra_pre_name %>%
    left_join(cat_acc_deme, by=c("deme"="deme_of_origin"))
#str(cat_acc_merge_by_deme)

cat_acc_merge_by_deme %>%
    filter(is.na(latitude))

# Now there are just 8 missing records, for the two demes with no lat/long
# info.  Let's use this set going forward

cat_acc_merge_by_deme_pre = cat_acc_merge_by_deme %>%
    select(sample_name, locality, latitude, longitude) %>%
    mutate(country="Spain") %>%
    filter(!is.na(locality))
#str(cat_acc_merge_by_deme_pre)

cat_meta = cat_sra_pre %>%
    inner_join(cat_acc_merge_by_deme_pre, by="sample_name") %>%
    mutate(oa_id=sprintf("OACT%04d", 1:n()))
str(cat_meta)

all_meta = bind_rows(all_meta, cat_meta)


# ## Tibetian
#
# [Zeng et al. 2017](https://doi.org/10.1016/j.scib.2017.10.007)
#
# No metadata table, but in the sups they say:
#
# > The seeds and leaves of Tibet-0 were collected in the wild forest of
# > Duilongdeqing County (N29.6903, E90.9338, altitude: 4200m asl), Tibet in 2013
# > when it had blossomed and borne fruit. 

tib_sra = read_csv("source/tibetian/SraRunTable.txt")
#str(tib_sra[])
tib_sra_pre = tib_sra %>%
    select(sra_run=Run, bioproject=BioProject, sample_name=`Sample Name`,
           library_layout=LibraryLayout, instrument=Instrument,
           species=Organism)

tib_acc = tribble(
    ~sample_name, ~latitude, ~longitude, ~elevation, ~country,
    "Tibet_0", 29.6903, 90.9338, 4200, "Tibet",
)

tib_meta = full_join(tib_sra_pre, tib_acc, by="sample_name") %>%
    mutate(oa_id=sprintf("OATB%04d", 1:n()))
str(tib_meta)

all_meta = bind_rows(all_meta, tib_meta)


# ## Korean
# 
# [J. Park et al. 2020](https://www.hindawi.com/journals/ijg/2020/3236461/)
# One accession from Korea

kor_sra = read_tsv("source/korean/filereport_read_run_SRX10808800_tsv.txt")
#str(kor_sra[])
kor_sra_pre = kor_sra %>%
    select(sra_run=run_accession, bioproject=study_accession,
           sample_name=sample_alias, library_layout,
           instrument=instrument_model, species=scientific_name)

kor_acc = tribble(
    ~sample_name, ~latitude, ~longitude, ~collector, ~country, ~locality,
    "180404IB4", 35.242862, 126.508987, "Y. Kim", "South Korea", "Yeonggwang-gun",
)


kor_meta = full_join(kor_acc, kor_sra_pre, by="sample_name") %>%
    mutate(oa_id=sprintf("OAKR%04d", 1:n()))
str(kor_meta)
all_meta = bind_rows(all_meta, kor_meta)


# ## Novikova and friends
#
# A series of overlapping papers and BioProjects from mostly non-thaliana
# species of Arabidopsis. The highest profile paper was Polina Novikova's 2018
# Nat Genet paper, but there are several others from the Koch, Nordborg, and
# Shimizu labs. As these overlap I'm going to process these as a single block.

nov_acc_orig = read_xlsx("source/novikova/41588_2016_BFng3617_MOESM15_ESM.xlsx")
#str(nov_acc_orig)

nov_acc_pre = nov_acc_orig %>%
    transmute(
        species=clade,
        sample_name = sample, 
        latitude = lat, 
        longitude = long, 
        ploidy = ploidy,
        biosample = BioSample,
        bioproject = BioProject,
        oa_id = sprintf("OANK%04d", 1:n()),
    ) %>% 
    filter(!is.na(latitude), !is.na(longitude), species!= "outgroup") 

nov_bioproj = nov_acc_pre %>%
    pull(bioproject) %>%
    unique()

# Downloaded all metadata from ENA for the above bioprojs

nov_all_sra =  read_tsv("source/novikova/results_read_run_tsv.txt",
                        col_types=cols(collection_date=col_character()))
#str(nov_all_sra[])

nov_all_sra = nov_all_sra %>%
    mutate(matchable_sample_acc = ifelse(accession%in%nov_acc_pre$biosample,
                                         accession,
                                  ifelse(secondary_sample_accession%in%nov_acc_pre$biosample,
                                         secondary_sample_accession,
                                         NA)))


nov_sra_pre = nov_all_sra %>%
    filter(!is.na(matchable_sample_acc)) %>%
    select(sra_run=run_accession, sample_accession=matchable_sample_acc,
           library_layout, sample_alias, species=scientific_name,
           lat_sra = lat, lon_sra=lon, instrument_model,
           collector=collected_by, collection_date, 
)


nov_meta_pre  = left_join(nov_acc_pre, nov_sra_pre, by=c("biosample"="sample_accession")) 

# Some sanity checks: do the lats & longs from SRA and from the Novikova supps line up?
all(with(nov_meta_pre, latitude == lat_sra & longitude == lon_sra))
# then which don't?
nov_iffy = nov_meta_pre %>%
    filter(latitude != lat_sra | longitude != lon_sra)
nov_iffy


# Include only those samples whose lat/long match.
nov_meta_pre = nov_meta_pre %>%
    filter(is.na(lat_sra) | latitude == lat_sra,
           is.na(lon_sra) | longitude == lon_sra)
#str(nov_meta_pre)

nov_meta = nov_meta_pre %>%
    select(oa_id, sample_name, species=species.y, latitude, longitude, ploidy,
           collector, collection_date, sra_run, oa_id, biosample, bioproject,
           instrument_model, library_layout)

all_meta = bind_rows(all_meta, nov_meta)


# ## Paape et al.'s kamchatika
#
# Natural variation in A. kamchatika. From Paape et al. nat comms

pkm_acc_orig = read_tsv("source/paape-kamchatika/tabula-41467_2018_6108_MOESM1_ESM-just-accessions.tsv")
#str(pkm_acc_orig[])

pkm_acc_pre = pkm_acc_orig %>%
    mutate(latitude = parse_lat(latitude),
           longitude = parse_lon(longitude),
           oa_id = sprintf("OAPK%04d", seq_len(n())))

pkm_sra_orig = read_tsv("source/paape-kamchatika/filereport_read_run_PRJDB6166_tsv.txt")
#str(pkm_sra_orig[])

pkm_sra_pre = pkm_sra_orig %>%
    transmute(bioproject=study_accession, sra_run=run_accession,
              species=scientific_name, sample_name=library_name,
              library_layout, instrument_model, biosample=sample_accession)
        
pkm_meta_pre = full_join(pkm_acc_pre, pkm_sra_pre, by=c("biosample"))

# It turns out their metadata includes two accessions from Novikova et al above
# in a separate bioproject. Hence we have some accessions in the above missing
# SRA metadata. We want to exclude those, as we already have the samples
# gathered above. 

pkm_meta_pre = pkm_meta_pre %>%
    filter(!is.na(bioproject))

# there are also some duplicated columns in the metadata, which we should check
# are equal and deduplicate.

all(with(pkm_meta_pre, sample_name.x == sample_name.y))
all(with(pkm_meta_pre, species.x == species.y))

pkm_meta = pkm_meta_pre %>%
    select(oa_id, sample_name=sample_name.x, latitude, longitude, country,
           locality=location, species=species.x, ploidy, sra_run, biosample,
           bioproject, library_layout, instrument_model)
str(pkm_meta[])

all(colnames(pkm_meta) %in% colnames(all_meta))

all_meta = bind_rows(all_meta, pkm_meta)



# ## Monnahan et al. A. arenosa
#
# from 10.1038/s41559-019-0807-4

mon_pops = read_csv("source/monnahan/tabula-41559_2019_807_MOESM1_ESM.csv") 
#str(mon_pops[])

mon_sra_orig = read_tsv("source/monnahan/filereport_read_run_PRJNA484107_tsv.txt")
#str(mon_sra_orig[])

mon_sra_pre = mon_sra_orig %>%
    select(sra_run=run_accession, bioproject=study_accession,
           biosample=sample_accession, sample_name=sample_alias,
           library_layout, instrument=instrument_model,
           species=scientific_name)

# So the "accession" metadata is acutally by sampling locality, so we need to
# extract the population codes out of the sample names then match on that.

mon_meta_pre = mon_sra_pre %>%
    mutate(Pop = substr(sample_name, 1,3)) %>%
    left_join(mon_pops, by="Pop") %>%
    mutate(observed_nind = table(Pop)[Pop])

# All SRAs have pop-level metadata?
all(!is.na(mon_meta_pre$latitude))

# All pops have the correct N indiv?
all(with(mon_meta_pre, observed_nind == Nind))

# The country codes are all non-standard, we should fix them.
table(mon_meta_pre$Country)
country_fix = c("Sk"="Slovakia", "Pol"="Poland", "At"="Austria",
                "BH"="Bosnia and Herzegovina", "Srb"="Serbia",
                "Cro"="Croatia", "Bel"="Belgium", "Rom"="Romania",
                "Hun"="Hungary", "Swe"="Sweden", "CZ"="Czech Republic",
                "D"="Germany", "Lit"="Lithuania")

mon_meta = mon_meta_pre %>%
    mutate(oa_id=sprintf("OAAN%04d", seq_len(n())),
           country=as.character(country_fix[Country])) %>%
    separate(Locality, c("locality", "collection_notes"), sep=", ", remove=T,
             extra="merge", fill="right") %>%
    select(oa_id, sample_name, species, latitude, longitude, elevation, country,
           locality, collection_notes, sra_run, bioproject, biosample, library_layout,
           instrument, ploidy)

all_meta = bind_rows(all_meta, mon_meta)

# ## Lucek & Willi North American A. lyrata
#
# From 10.1371/journal.pgen.1009477
#
# This metadata doesn't seem complete, and by population. However, in the SRA's
# run exporter, all relevant metadata shows up, and seems correct and complete.
# So, we use the metadata from the SRA, but do a quick check of lats and longs
# of the data that is there from the paper to check it all concurrs.

lwl_pop = read_xlsx("source/lucek_willi/journal.pgen.1009477.s008.xlsx", skip=2)
#str(lwl_pop)

lwl_sra = read_csv("source/lucek_willi/SraRunTable.txt")
#str(lwl_sra[])

lwl_meta = lwl_sra %>%
    select(sra_run=Run, bioproject=BioProject, biosample=BioSample,
           country=geo_loc_name_country, instrument_model=Instrument,
           lat_lon, library_layout=LibraryLayout, species=Organism, sample_name) %>%
    mutate(oa_id=sprintf("OALL%04d", as.numeric(as.factor(biosample)))) %>%
    select(oa_id, sample_name, everything()) %>%
    tidyr::extract(lat_lon, into=c("latitude", "longitude"),
                   regex="(\\d+\\.\\d+ ?[NS]) ?(\\d+\\.\\d+ ?[EW])") %>%
    mutate(latitude = parse_lat(latitude),
           longitude = parse_lon(longitude))
#str(lwl_meta)

lwl_pop = lwl_pop %>%
    transmute(pop = Population, latitude_supps=`Latitude (°N)`,
           longitude_supps=-`Longitude (°W)`)

lwl_check = lwl_meta %>%
    mutate(pop = str_split(sample_name, "_", n=4, simplify=T)[,3]) %>%
    left_join(lwl_pop, by=c("pop"))

lwl_bad = lwl_check %>%
    filter(!is.na(latitude_supps) & round(latitude_supps, 2) != latitude |
           !is.na(longitude_supps) & round(longitude_supps, 2) != longitude)
str(lwl_bad)

# So in short, the lat/longs from the supps are a lot more precise, but when
# rounded down to two digits (like the ones in the SRA), they match exactly for
# all records.
#
# Since we went to the effort of extracting the data from the supps and
# matching it up, I'm going to use the lat/longs from the supps unless they're
# not present, in which case we'll use the SRA ones.

#str(lwl_check)
lwl_meta_pop = lwl_check %>%
    mutate(latitude = ifelse(!is.na(latitude_supps), latitude_supps, latitude),
           longitude = ifelse(!is.na(longitude_supps), longitude_supps, longitude)) %>%
    select(-pop, -latitude_supps, -longitude_supps)

str(lwl_meta_pop)

all_meta = bind_rows(all_meta, lwl_meta_pop)

# ## Lee et al. A. halleri
#
# From a biorxiv paper, 10.1101/859249
# The Wall_07 and Wall_10 samples seem missing from SRA?

lh_sra = read_tsv("source/halleri-lee-etal/filereport_read_run_PRJEB35573_tsv.txt")
#str(lh_sra[])

lh_sra_pre = lh_sra %>%
    filter(library_strategy != "RNA-Seq") %>%
    select(sra_run=run_accession, bioproject=study_accession,
           biosample=sample_accession, sample_name=sample_alias,
           library_layout, instrument=instrument_model,
           species=scientific_name)
#str(lh_sra_pre)

lh_acc = read_csv("source/halleri-lee-etal/tabula-fixed-supps-accessions.csv")

lh_acc_pre = lh_acc %>%
    mutate(sample_name = gsub("_0?", "", sample_name))

lh_meta = inner_join(lh_sra_pre, lh_acc_pre, by="sample_name") %>%
    mutate(oa_id=sprintf("OALH%04d", as.numeric(as.factor(biosample)))) %>%
    mutate(latitude = parse_lat(latitude), longitude = parse_lon(longitude))
str(lh_meta)

all_meta = bind_rows(all_meta, lh_meta)


# ## Frachon et al. Pyrnees
#
# This is poolseq unfortunately, but might be useful anyway as it's very nice
# dense sampling in a moderately sized region.

fp_sites = read_csv("source/frachon-midi-pyrnees/all-meta-imgextracted.csv") %>%
    janitor::clean_names()
#str(fp_sites[])

fp_sra = read_csv("source/frachon-midi-pyrnees/SraRunTable.txt")
#str(fp_sra[])

fp_sra_pre = fp_sra %>%
    select(sra_run=Run, bioproject=`SRA Study`, biosample=BioSample,
           instrument_model=Instrument, library_layout=LibraryLayout,
           species=Organism, sample_name=title) %>%
    mutate(sample_name = gsub("^Climares_Ath_", "", sample_name),
           # There is a second pool for the MONTM-A population, MONTM-Abis. So
           # we regex a way the bis suffix so it matches.
           sample_name_match = gsub("bis$", "", sample_name))
#str(fp_sra_pre)

# Do all sites match the lat/longs?
all(with(fp_sra_pre, sample_name_match %in% fp_sites$population))

fp_meta = left_join(fp_sra_pre, fp_sites, by=c("sample_name_match"="population")) %>%
    mutate(oa_id=sprintf("OAFP%04d", as.numeric(as.factor(biosample)))) %>%
    select(oa_id, sample_name, latitude, longitude, elevation=altitude,
           locality=town, species, sra_run, bioproject, biosample,
           library_layout, instrument_model, pooled_plants=no_of_pooled_plants)
str(fp_meta)

all_meta = bind_rows(all_meta, fp_meta)


# ## Fulgione et al.  Maiderans
#
# from 10.1093/molbev/msx300
#
# The data isn't on the sra, but we have the reads internally

fm_acc = read_csv("source/madeirans/msx300_supp.csv")
#str(fm_acc[])

fm_meta = fm_acc %>%
    mutate(oa_id=sprintf("OAFM%04d", as.numeric(as.factor(ID))), internal_data="Y") %>%
    select(oa_id, sample_name=ID, locality=Region, latitude=Latitude,
           longitude=Longitude, elevation=Altitude, collector=Collector,
           internal_data)
str(fm_meta)

all_meta = bind_rows(all_meta, fm_meta)

# ## TODO Consistent unique IDs for all accessions
#
# Not all plants have ecotype IDs, and many names have charachers we dont' want
# in file names (non-ascii, or weird punctuation). For our general sanity, we
# want a simple ID like an ecotype ID for all samples. For samples with ecotype
# IDs, we'll just use the ecotype ID. For others, we'll make new ones.

str(all_meta)

all_sra = all_meta %>%
    select(sra_run, bioproject, instrument, library_layout, oa_id) %>%
    unique()

all_acc = all_meta %>%
    select(oa_id, ecotype_id, sample_name, species, latitude, longitude,
           elevation, locality, country, collector, collection_date, cs_number) %>%
    mutate(dataset=substr(oa_id, 3, 4)) %>%
    unique()


# Did we miss any columns?
setdiff(union(colnames(all_acc), colnames(all_sra)), colnames(all_meta))

str(all_acc)
str(all_sra)

table(all_acc$dataset)


# # Validate metadata


# ## Sra metadata validation
#
# This is pretty cursary, chekcing that the run IDs look legit, and there isn't
# any duplication.

sra_val = validator(
    run = grepl("(ERR|SRR)\\d+", sra_run, perl=T),
    run_unq = multiplicity(sra_run) == 1,
    bioproj = !is.na(bioproject),
    liblay = library_layout %in% c("PAIRED", "SINGLE")
)

all_sra_val = confront(all_sra, sra_val)
summary(all_sra_val)

all_sra_bad = violating(all_sra, all_sra_val)
all_sra_bad

# ## Accession metadata

acc_val = validator(
    oa_unq = multiplicity(oa_id) == 1,
    eco_unq = multiplicity(ecotype_id) == 1,
    # Some names are duplicated, but have unique ecotype IDs.
    name_unq  = (!is.na(ecotype_id)) | multiplicity(sample_name) == 1,
    latlong_ok = !(is.na(latitude) | is.na(longitude)),
    dist2land = dist_to_land_km(longitude, latitude) < 5,
    spp = !is.na(species) | !grepl("^Arabidopsis", species)
)

all_acc_val = confront(all_acc, acc_val)
summary(all_acc_val)
errors(all_acc_val)

all_acc_bad = all_acc %>%
    mutate(failed = get_failing_tests(all_acc_val)) %>%
    filter(failed != "")

#write_tsv(all_acc_bad, "all_acc_bad.tsv")


# ## Corrections to metadata as a result of validation
#
# Some issues have been higlighted by the validation process. The following
# code fixes these in a documented way.

# One afican sample has the latitude up-side down. It's supposedly in south
# africa, but appears int the middle of the medeterrainian sea. Inverting the
# latitude 

all_acc =  all_acc %>%
    rows_update(tribble(
        ~oa_id, ~latitude, ~longitude,
        "OAAF0072", -33.399, 19.282
        ), by="oa_id")


# ## Final metadata validation
#
# After the fixes we applied above, is everything kosher?

all_acc_val = confront(all_acc, acc_val)
summary(all_acc_val)
errors(all_acc_val)

all_acc_bad = all_acc %>%
    mutate(failed = get_failing_tests(all_acc_val)) %>%
    filter(failed != "")

write_tsv(all_acc_bad, "all_acc_still_bad.tsv")
write_tsv(all_acc, "omniath_all_accessions.tsv")
write_tsv(all_sra, "omniath_all_sra.tsv")


# # Mapping of accessions

globe_bbox = c(left=-170, bottom=-58, right=179.99, top=75)
baselayer = get_stamenmap(globe_bbox, zoom = 3, maptype = "toner-background")

ggmap(baselayer) +
    geom_point(aes(x=longitude, y=latitude, colour=dataset), data=all_acc)


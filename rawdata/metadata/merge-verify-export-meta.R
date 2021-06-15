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

library(CoordinateCleaner)
library(countrycode)
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
#                   "readxl", "janitor", "countrycode", "CoordinateCleaner"))


match_country = function(eastings, northings, ret, crs=4326) {
    world = sf::st_as_sf(rworldmap::getMap(resolution="high"))
    everything = data.frame(eastings, northings, i=seq_along(eastings))
    points = everything %>%
        filter(!is.na(eastings), !is.na(northings)) %>%
        sf::st_as_sf(coords=1:2, crs=crs)
    points$nearest = sf::st_nearest_feature(points, world)
    points$ISO3 = world[points$nearest,]$ISO3
    points$dist = sf::st_distance(points, world[points$nearest,], by_element=T) %>%
        units::set_units("km") %>%
        as.numeric()
    everything = dplyr::left_join(everything, points, by='i')
    return(everything[[ret]])
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

#sra.sql = "/home/kevin/data/work/2021-06-06_SRAmetadb.sqlite.gz"
#sra = dbConnect(SQLite(), sra.sql)

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
    select(ecotype_id=Ecotype, sra_run=Run, biosample=BioSample,
           bioproject=BioProject, library_layout=LibraryLayout,
           instrument_model=Instrument, species=Organism)
stopifnot(all(kg_sra_pre$ecotype_id %in% kg_acc_pre$ecotype_id))
stopifnot(all(kg_acc_pre$ecotype_id %in% kg_sra_pre$ecotype_id))
#str(kg_sra_pre)


kg_meta = inner_join(kg_acc_pre, kg_sra_pre, by="ecotype_id")
str(kg_meta)

all_meta = kg_meta


# ## Cao et al Resequencing
#
# This is fernando's re-sequecning of the lines from Cao et al. 2011.

cr_acc = read_tsv("source/internal-data/cao-reseq.txt")

cr_meta = kg_acc_pre %>%
    select(-oa_id) %>%
    inner_join(cr_acc, by="ecotype_id") %>%
    mutate(oa_id = sprintf("OACR%04d", ecotype_id),
           internal_data=T, internal_sample_name=as.character(ecotype_id))

all_meta = bind_rows(all_meta, cr_meta)


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
    select(sra_run=Run, bioproject=BioProject, biosample=BioSample,
           sample_name=Alias, library_layout=LibraryLayout,
           instrument_model=Instrument, species=Organism)
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

setdiff(colnames(af_meta_complete), colnames(all_meta))
all_meta = bind_rows(all_meta, af_meta_complete)


# ## Argentinian lines
#
# From [Luciana Kasulin's paper in MolEcol](https://dx.doi.org/10.1111/mec.14107)

arg_sra = read_tsv("source/argentinians/SraRunTable_PRJEB9862.txt")
#str(arg_sra[])

arg_sra_pre = arg_sra %>%
    select(sra_run=run_accession, bioproject=study_accession,
           biosample=sample_accession, sample_name=sample_alias,
           library_layout, instrument_model, species=scientific_name)
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

setdiff(colnames(arg_meta), colnames(all_meta))
all_meta = bind_rows(all_meta, arg_meta)


# ## Chinese lines
#
# From Yangtze river basin, see [Zou et al.
# 2017](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1378-9).

chi_sra = read_tsv("source/chinese/SraRunTable_PRJNA293798.txt")
#str(chi_sra[])
chi_sra_pre = chi_sra %>%
    select(sra_run=run_accession, bioproject=study_accession,
           biosample=sample_accession, sample_name=sample_alias,
           library_layout, instrument_model, species=scientific_name)
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
setdiff(colnames(chi_meta), colnames(all_meta))
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
    select(sra_run=Run, bioproject=BioProject, biosample=BioSample,
           sample_name=Sample_Name, library_layout=LibraryLayout,
           instrument_model=Instrument, species=Organism)
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

setdiff(colnames(cat_meta), colnames(all_meta))
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
    select(sra_run=Run, bioproject=BioProject, biosample=BioSample,
           sample_name=`Sample Name`, library_layout=LibraryLayout,
           instrument_model=Instrument, species=Organism)

tib_acc = tribble(
    ~sample_name, ~latitude, ~longitude, ~elevation, ~country,
    "Tibet_0", 29.6903, 90.9338, 4200, "Tibet",
)

tib_meta = full_join(tib_sra_pre, tib_acc, by="sample_name") %>%
    mutate(oa_id=sprintf("OATB%04d", 1:n()))
str(tib_meta)

setdiff(colnames(tib_meta), colnames(all_meta))
all_meta = bind_rows(all_meta, tib_meta)


# ## Korean
# 
# [J. Park et al. 2020](https://www.hindawi.com/journals/ijg/2020/3236461/)
# One accession from Korea

kor_sra = read_tsv("source/korean/filereport_read_run_SRX10808800_tsv.txt")
#str(kor_sra[])
kor_sra_pre = kor_sra %>%
    select(sra_run=run_accession, bioproject=study_accession,
           biosample=sample_accession, sample_name=sample_alias,
           library_layout, instrument_model, species=scientific_name)

kor_acc = tribble(
    ~sample_name, ~latitude, ~longitude, ~collector, ~country, ~locality,
    "180404IB4", 35.242862, 126.508987, "Y. Kim", "South Korea", "Yeonggwang-gun",
)


kor_meta = full_join(kor_acc, kor_sra_pre, by="sample_name") %>%
    mutate(oa_id=sprintf("OAKR%04d", 1:n()))
str(kor_meta)
setdiff(colnames(kor_meta), colnames(all_meta))
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
    ) %>% 
    unique() %>%
    mutate(oa_id = sprintf("OANK%04d", as.numeric(as.factor(biosample)))) %>%
    filter(!is.na(latitude), !is.na(longitude), species!= "outgroup") 

nov_bioproj = nov_acc_pre %>%
    pull(bioproject) %>%
    unique()

# Downloaded all metadata from ENA for the above bioprojs

nov_all_sra =  read_tsv("source/novikova/results_read_run_tsv.txt",
                        col_types=cols(collection_date=col_character()))
#str(nov_all_sra[])

nov_all_sra = nov_all_sra %>%
    mutate(matchable_sample_acc = case_when(
        accession %in% nov_acc_pre$biosample ~ accession,
        secondary_sample_accession %in% nov_acc_pre$biosample ~ secondary_sample_accession,
        TRUE ~ as.character(NA)
    ))

nov_sra_pre = nov_all_sra %>%
    filter(!is.na(matchable_sample_acc)) %>%
    select(
        sra_run=run_accession, sample_accession=matchable_sample_acc,
        library_layout, sample_alias, species=scientific_name, lat_sra = lat,
        lon_sra=lon, instrument_model, collector=collected_by, collection_date,
    )


nov_meta_pre  = left_join(nov_acc_pre, nov_sra_pre,
                          by=c("biosample"="sample_accession")) %>%
    unique() %>%
    mutate(mult = multiplicity(sra_run)) %>%
    arrange(mult, biosample)

# So the metadata in the supps table is actually per BAM, i.e. you get two
# records for each polyploid. We only want one record per sample, so we ablate
# the bits of the sample name that refer to the reference genome.

name_fix = c(
    "AkamchaticaKWShal"="AkamchaticaKWS",
    "AkamchaticaPAKhal"="AkamchaticaPAK",
    "AkamchaticaTOYhal"="AkamchaticaTOY",
    "AkamchaticaTWNhal"="AkamchaticaTWN",
    "AkamchaticaKWSlyr"="AkamchaticaKWS",
    "AkamchaticaPAKlyr"="AkamchaticaPAK",
    "AkamchaticaTOYlyr"="AkamchaticaTOY",
    "AkamchaticaTWNlyr"="AkamchaticaTWN",
    "Aar.AS459"="AS459",
    "Aar.ASO5" ="ASO5",
    "Aar.ASS3a"="ASS3a",
    "Ath.AS459"="AS459",
    "Ath.ASO5" ="ASO5",
    "Ath.ASS3a"="ASS3a"
)

nov_meta_pre = nov_meta_pre %>%
    mutate(sample_name = ifelse(sample_name %in% names(name_fix),
                                name_fix[sample_name],
                                sample_name)) %>%
    select(-mult) %>%
    unique() %>%
    mutate(mult = multiplicity(sra_run))


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

setdiff(colnames(nov_meta), colnames(all_meta))
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

setdiff(colnames(pkm_meta), colnames(all_meta))
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
           library_layout, instrument_model,
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
           instrument_model, ploidy)

setdiff(colnames(mon_meta), colnames(all_meta))
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

setdiff(colnames(lwl_meta_pop), colnames(all_meta))
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
           library_layout, instrument_model,
           species=scientific_name)
#str(lh_sra_pre)

lh_acc = read_csv("source/halleri-lee-etal/tabula-fixed-supps-accessions.csv")

lh_acc_pre = lh_acc %>%
    mutate(sample_name = gsub("_0?", "", sample_name)) %>%
    select(-pop_code, -pop_size)

lh_meta = inner_join(lh_sra_pre, lh_acc_pre, by="sample_name") %>%
    mutate(oa_id=sprintf("OALH%04d", as.numeric(as.factor(biosample)))) %>%
    mutate(latitude = parse_lat(latitude), longitude = parse_lon(longitude))
str(lh_meta)

setdiff(colnames(lh_meta), colnames(all_meta))
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
    mutate(oa_id=sprintf("OAFP%04d", as.numeric(as.factor(biosample))),
           country="France") %>%
    select(oa_id, sample_name, latitude, longitude, elevation=altitude,
           locality=town, species, sra_run, bioproject, biosample,
           library_layout, instrument_model, pooled_plants=no_of_pooled_plants,
           country)
str(fp_meta)

setdiff(colnames(fp_meta), colnames(all_meta))
all_meta = bind_rows(all_meta, fp_meta)


# ## Fulgione et al.  Maiderans
#
# from 10.1093/molbev/msx300
#
# The data isn't on the sra, but we have the reads internally

fm_acc = read_csv("source/madeirans/msx300_supp.csv")
#str(fm_acc[])

fm_meta = fm_acc %>%
    mutate(oa_id=sprintf("OAFM%04d", as.numeric(as.factor(ID))), internal_data=T,
           country="Spain") %>%
    select(oa_id, sample_name=ID, locality=Region, latitude=Latitude,
           longitude=Longitude, elevation=Altitude, collector=Collector,
           internal_data, internal_sample_name=ID, country)
str(fm_meta)

setdiff(colnames(fm_meta), colnames(all_meta))
all_meta = bind_rows(all_meta, fm_meta)


# ## Koch lab 
#
# Two different papers and bioprojects on lyrata and arenosa hybrid zones
#
# The first, Hohmann & Koch, is kl
# https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-017-4220-6

hl_sra = read_csv("source/koch_lab_lyrata/SraRunTable_PRJEB34247.txt",
                   col_types=cols(collection_date=col_character())) %>%
    janitor::clean_names()

hl_meta = hl_sra %>%
    mutate(oa_id=sprintf("OAHL%04d", as.numeric(as.factor(sample_name)))) %>%
    select(oa_id, sra_run=run, sample_name, bioproject=bio_project,
           biosample=bio_sample, collector=collected_by, collection_date,
           country=geo_loc_name_country, latitude=geographic_location_latitude,
           longitude=geographic_location_longitude,
           locality=geographic_location_region_and_locality,
           instrument_model=instrument, library_layout, species=organism)

# The second one is from Marburger et al. Nat comms.

ml_sra = read_csv("source/koch_lab_lyrata/SraRunTable_PRJEB20573.txt")
#str(ml_sra[])

ml_sra_pre = ml_sra %>%
    select(sra_run=Run, sample_name=Alias, bioproject=BioProject,
           biosample=BioSample, instrument_model=Instrument,
           library_layout=LibraryLayout, species=Organism) %>%
    # In the sra all the sample names have a " WGS" suffix, remove this
    mutate(sample_name =sub(" WGS$", "", sample_name))
#str(ml_sra_pre)

ml_acc = read_xlsx("source/koch_lab_lyrata/12864_2017_4220_MOESM2_ESM.xlsx", 
                   skip=1) %>%
    janitor::clean_names()

ml_acc_pre =  ml_acc %>%
    unite("species_parsed", genus, species, subspecies, sep=" ", na.rm=T, remove=F) %>%
    mutate(ploidy=as.character(ploidy)) %>%
    filter(bio_project=="PRJEB20573") %>%
    mutate(country=str_split(locality, ";", simplify=T)[,1],
           herbarium_id=sample_id_herbarium_code_heid_heidelberg_uzh_zurich,
           sample_name=sub("^[^:].+::", "", herbarium_id),
           locality = sub("\\n", " ", locality),
           latitude= parse_lat(gsub(",", ".", latitude)),
           longitude= parse_lon(gsub(",", ".", longitude)),
           oa_id=sprintf("OAML%04d", as.numeric(as.factor(sample_name)))) %>%
    select(oa_id, sample_name, species=species_parsed, ploidy, locality,
           latitude, longitude, collector, collection_date,
           bioproject=bio_project, country)
str(ml_acc_pre)

ml_meta = full_join(ml_acc_pre, ml_sra_pre)

setdiff(colnames(ml_meta), colnames(all_meta))
setdiff(colnames(hl_meta), colnames(all_meta))
all_meta = bind_rows(all_meta, ml_meta, hl_meta)


# ## Serpentine soil halleri and arenosa
#
# From [a paper from Levi Yant's
# lab](https://royalsocietypublishing.org/doi/full/10.1098/rstb.2018.0243)

sh_pops = read_tsv("source/arenosa-halleri-PRJNA506705/pops.tsv") %>%
    clean_names() %>%
    mutate(latitude=parse_lat(latitude), longitude=parse_lon(longitude))
#str(sh_pops[])

sh_sra = read_csv("source/arenosa-halleri-PRJNA506705/SraRunTable.txt")
#str(sh_sra[])

sh_sra_pre = sh_sra %>%
    select(sra_run=Run,  bioproject=BioProject, biosample=BioSample,
           country=geo_loc_name_country, instrument_model=Instrument,
           sample_name=`Library Name`, library_layout=LibraryLayout,
           species=Organism, lat_lon) %>%
    mutate(popcode=str_split(sample_name, "_", n=2, simplify=T)[,1]) %>%
    mutate(oa_id=sprintf("OASH%04d", as.numeric(as.factor(sample_name)))) %>%
    arrange(popcode, oa_id)
#str(sh_sra_pre)

sh_meta = left_join(sh_sra_pre, sh_pops, by=c("popcode"="population"))  %>%
    select(-lat_lon, -popcode, -soil_type) %>%
    rename(elevation=altitude)

str(sh_meta)

setdiff(colnames(sh_meta), colnames(all_meta))
all_meta = bind_rows(all_meta, sh_meta)


# ## PRJEB23202 
#
# This seems like a nice dataset, but I can't find the paper. Submission
# information suggests it's from Yvonne Willi's lab at Basel.

wn_sra = read_csv("source/PRJEB23202/SraRunTable.txt") %>%
    clean_names()
#str(wn_sra[])

wn_meta_pre = wn_sra %>%
    select(sra_run=run, sample_name=alias, bioproject=bio_project,
           biosample=bio_sample, collection_date, country=geo_loc_name_country,
           lat_lon, instrument_model=instrument, library_layout,
           species=organism) %>%
    mutate(
        oa_id = sprintf("OAWN%04d", as.numeric(as.factor(sample_name))),
        lat_lon = sub("(\\d+\\.\\d+)\\s*([NSEW])\\s*(\\d+\\.\\d+)\\s*([NSEW])",
                      "\\2\\1 \\4\\3", lat_lon, perl=T),
    )

# This is a bit of a faff, but it's how the parzer package extracts lat&long
wn_ll = parse_llstr(wn_meta_pre$lat_lon) %>% as.matrix()
wn_ll[!is.finite(wn_ll)] = NA
wn_ll = as.data.frame(wn_ll)

wn_meta = wn_meta_pre %>%
    mutate(longitude = wn_ll$lon, latitude=wn_ll$lat) %>%
    select(-lat_lon)
str(wn_meta)

# Any extra cols?
setdiff(colnames(wn_meta), colnames(all_meta))
all_meta = bind_rows(all_meta, wn_meta)

        
# ## Polish arenosa
#
# PRJNA667586, from [Konečná et
# al.](https://www.biorxiv.org/content/10.1101/2021.01.15.426785v1)

ka_sra = read_csv("source/arenosa-PRJNA667586/SraRunTable.txt") %>%
    clean_names()
#str(ka_sra[])

ka_pops = read_csv("source/arenosa-PRJNA667586/media-1-2.pdf.csv") %>%
    clean_names()
#str(ka_pops[])

ka_pops_pre = ka_pops %>%
    select(soil_type=bedrock, latitude=lat, longitude=lon, ploidy,
           elevation=altitude, pop, pop_code, locality=pop_name) %>%
    mutate(collection_notes = sprintf("soil_type=%s", soil_type))

ka_sra_pre = ka_sra %>%
    select(sample_name, population, biosample=bio_sample,
           bioproject=bio_project, sra_run=run, library_layout,
           instrument_model=instrument, country=geo_loc_name_country,
           locality=geo_loc_name, species=organism) %>%
    mutate(locality=sub("^[^:]+: ", "", locality))
#str(ka_sra_pre)

ka_meta_pre = ka_sra_pre %>%
    left_join(ka_pops_pre, by=c("population"="pop", "locality")) %>%
    filter(!is.na(latitude)) %>%
    select(-population, -pop_code, -soil_type) %>%
    mutate(oa_id = sprintf("OAKA%04d", as.numeric(as.factor(biosample))))
#str(ka_meta_pre)

setdiff(colnames(ka_meta_pre), colnames(all_meta))

all_meta = bind_rows(all_meta, ka_meta_pre)


# ## Hamala et al lyrata
#
# [Hamala & Savolainen](https://doi.org/10.1093/molbev/msz149) is the main
# paper, but some lats and longs are only in a different paper from the same
# lab, [Pyhäjärvi et al.](https://doi.org/10.3732/ajb.1100580).

hm_sra_1 = read_csv("source/hamala-lyrata/SraRunTable.txt") %>%
    clean_names()
hm_sra_2 = read_csv("source/hamala-lyrata/SraRunTable_PRJNA459481.txt") %>%
    clean_names()

hm_sra = bind_rows(hm_sra_1, hm_sra_2)
#str(hm_sra[])

hm_pop = read_csv("source/hamala-lyrata/page-1_table-1.csv") %>%
    mutate(latitude=parse_lat(latitude), longitude=parse_lon(longitude)) %>%
    mutate(geo_loc_name = sprintf("%s: %s", country, locality))
#str(hm_pop)

hm_meta_pre = left_join(hm_sra, hm_pop, by="geo_loc_name")
#str(hm_meta_pre[])

hm_meta = hm_meta_pre %>%
    select(sra_run=run, bioproject=bio_project, biosample=bio_sample,
           instrument_model=instrument, sample_name, species=organism, 
           country, locality, latitude, longitude, collection_notes, collector) %>%
    mutate(oa_id = sprintf("OAHS%04d", as.numeric(as.factor(biosample))))

setdiff(colnames(hm_meta), colnames(all_meta))
all_meta = bind_rows(all_meta, hm_meta)


# ## Parker et al. Snowdonia
#
# [Parker et al.](https://dx.doi.org/10.1038/s41598-017-08461-5)

ps_sra = read_csv("source/snowdonia/SraRunTable.txt",
                  col_types=cols(collection_date=col_character())) %>%
    clean_names()
#str(ps_sra[])

ps_meta_pre = ps_sra %>%
    select(sra_run=run, sample_name=title, bioproject=bio_project,
           biosample=bio_sample, collector=collected_by, collection_date,
           country=geo_loc_name_country,
           elevation=geographic_location_altitude,
           latitude=geographic_location_latitude,
           longitude=geographic_location_longitude,
           locality=geographic_location_region_and_locality,
           instrument_model=instrument, library_layout, species=organism,
           ploidy) %>%
    mutate(oa_id = sprintf("OAPS%04d", as.numeric(as.factor(biosample)))) %>%
    filter(instrument_model != "MinION")

setdiff(colnames(ps_meta_pre), colnames(all_meta))
str(ps_meta_pre)

all_meta = bind_rows(all_meta, ps_meta_pre)


# # Guggisberg Bohemian Arabidopsis
#
# A. thaliana, negelecta, and lyrata from Bohemia

bo_sra = read_csv("source/guggisberg/SraRunTable.txt",
                  col_types=cols(collection_date=col_character())) %>%
    clean_names()
bo_acc = read_tsv("source/guggisberg/mec14930-sup-0002-tables1.csv",
                  col_types=cols(Date=col_character())) %>%
    clean_names()

bo_sra_pre = bo_sra %>%
    transmute(
        sra_run=run, sample_name=isolate, bioproject=bio_project,
        biosample=bio_sample, collector=collected_by, collection_date,
        collection_notes=sprintf("soil_type=%s", isolation_source),
        country=geo_loc_name_country, locality=geo_loc_name,
        instrument_model=instrument, library_layout, species=organism,
        lat_lon
    ) %>%
    mutate(oa_id = sprintf("OABO%04d", as.numeric(as.factor(biosample)))) %>%
    mutate(lat_lon = sub("(N|S) ", "\\1, ", lat_lon) %>%
                     parse_llstr(),
           lat_sra = lat_lon[,1],
           lon_sra = lat_lon[,2],
    ) %>% select(-lat_lon)

bo_meta = full_join(bo_sra_pre, bo_acc, by=c("sample_name"="acronym")) %>%
    mutate(latitude = ifelse(is.na(latitude), lat_sra, latitude),
           longitude = ifelse(is.na(longitude), lon_sra, longitude)) %>%
    select(oa_id, sra_run, sample_name, bioproject, biosample, collector,
           collection_date, country=country.x, locality, instrument_model,
           latitude, longitude, library_layout, elevation=altitude,
           species=species.x)

setdiff(colnames(bo_meta), colnames(all_meta))
str(bo_meta)

all_meta = bind_rows(all_meta, bo_meta)


# ## Gould et al North american Ath
#
# PRJNA288374

gn_sra = read_csv("source/gould-north-american/SraRunTable_PRJNA288374.txt") %>%
    clean_names()
#str(gn_sra[])

gn_acc = read_tsv("source/gould-north-american/mec13643-sup-0002-tables1.csv",
                  na="--") %>% clean_names()
#str(gn_acc[])

gn_sra_pre = gn_sra %>%
    select(sra_run=run, sample_name, bioproject=bio_project,
           biosample=bio_sample, country=geo_loc_name_country,
           instrument_model=instrument, library_layout, species=organism)

gn_acc_pre = gn_acc %>%
    select(sample_name=line_name, latitude, longitude, cs_number) %>%
    mutate(longitude = -longitude) # these are degrees W

gn_meta = left_join(gn_sra_pre, gn_acc_pre) %>%
    mutate(oa_id = sprintf("OAGN%04d", as.numeric(as.factor(biosample)))) 

setdiff(colnames(gn_meta), colnames(all_meta))
all_meta = bind_rows(all_meta, gn_meta)


# # Shirsekar et al North American lines
#
# So this is pretty incomplete, but I'll follow it up with Gautam once his
# revisions are done.

sn_sra = read_csv("source/shirsekar-north-american/SraRunTable.txt") %>%
    clean_names() %>%
    filter(library_selection != "Reduced Representation")
#str(sn_sra[])

sn_sra_pre = sn_sra %>%
    transmute(sra_run=run, sample_name_sra=sample_name_2, bioproject=bio_project,
              biosample=bio_sample, collection_date=as.character(collection_year),
              instrument_model=instrument, library_layout, species=organism) %>%
    mutate(sample_name = sub("^r\\d+_l\\d+_", "", sample_name_sra))
#str(sn_sra_pre[])

sn_acc = read_tsv("source/shirsekar-north-american/srisekar-acc.csv") %>% 
    clean_names()
#str(sn_acc[])

sn_acc_pre = sn_acc %>%
    transmute(sample_name = sub("^ind_", "", individual), latitude, longitude,
              country, locality=population)
#str(sn_acc_pre[])

sn_meta = inner_join(sn_sra_pre, sn_acc_pre, by="sample_name") %>%
    mutate(oa_id = sprintf("OASN%04d", as.numeric(as.factor(biosample))))  %>%
    select(-sample_name_sra)
str(sn_meta)

setdiff(colnames(sn_meta), colnames(all_meta))
all_meta = bind_rows(all_meta, sn_meta)


# ## Günther et al Tyrolian
#
# 5 poolseq pops from alpine and valley pops in Tyrol

gt_sra = read_csv("source/günther-alpine-italy/SraRunTable_PRJEB12316.txt") %>%
    clean_names()
#str(gt_sra[])

gt_sra_pre = gt_sra %>%
    select(sra_run=run, sample_name=alias, bioproject=bio_project,
           biosample=bio_sample, instrument_model=instrument, library_layout,
           species=organism)

gt_pops = read_tsv("source/günther-alpine-italy/sites.tsv") %>%
    clean_names()
#str(gt_pops[])

gt_pops_pre = gt_pops %>%
    transmute(population=sub("/", "", population), pooled_plants=individuals,
              elevation, latitude=parse_lat(latitude),
              longitude=parse_lon(longitude), country="ITA",
              oa_id = sprintf("OAGT%04d", as.numeric(as.factor(population))))


gt_meta_pre = inner_join(gt_sra_pre, gt_pops_pre, by=c("sample_name"="population"))
str(gt_meta_pre)

setdiff(colnames(gt_meta_pre), colnames(all_meta))
all_meta = bind_rows(all_meta, gt_meta_pre)


# # Global dataset modification
# 
# There are some modifications that need to be done for the whole dataset, now
# that we have imported all individual datasets. There are also some fixes here
# for problems previously identified in the validation process below.

# ## Normalise country codes
#
# There are a whole panoply of ways of refering to countries in the `country`
# column of our metadata. The below code normalises this to two columns, an ISO
# 3-letter code, and a normalised country name.

custom_codes = c(
    "GER"="DEU",
    "NED"="NLD",
    "SUI"="CHE",
    "POR"="PRT",
    "DEN"="DNK",
    "UNK"=NA,
    "CRO"="HRV",
    "BUL"="BGR",
    "Tibet"="CHN"
)

ctry = tibble(orig=na.omit(unique(all_meta$country))) %>%
    mutate(
        country_code = case_when(
            (orig %in% names(custom_codes)) ~ custom_codes[orig],
            (orig %in% codelist$iso2c) ~ countrycode(orig, "iso2c", "iso3c"),
            (orig %in% codelist$iso3c) ~ orig,
            TRUE ~ countrycode(orig, "country.name.en", "iso3c")
        ),
        country_name=countrycode(country_code, "iso3c", "country.name"))

all_meta = all_meta %>%
    left_join(ctry, by=c("country"="orig"))

# ## Normalise ploidy values
#
# Likewise, ploidies are refered to in a number of different codes. We sanitise
# all to words, e.g. diploid, tetraploid, allotetraploid, etc. 4x is asssumed
# to be autotetraploid.

all_meta = all_meta %>%
    mutate(ploidy = ifelse(grepl("2", ploidy), "diploid",
                    ifelse(grepl("4", ploidy), "tetraploid",
                    ifelse(grepl("allo", ploidy), "allotetraploid",
                           ploidy))))

# Which species don't have annotated ploidy?
all_meta %>%
    filter(is.na(ploidy)) %>%
    group_by(species) %>%
    summarise(n=n())

# ## Corrections to metadata as a result of validation
#
# Some issues have been higlighted by the validation process. The following
# code fixes these in a documented way.

# One afican sample has the latitude up-side down. It's supposedly in south
# africa, but appears in the middle of the medeterrainian sea. Inverting the
# latitude 

all_meta =  all_meta %>%
    mutate(latitude = case_when(
        oa_id=="OAAF0072" ~ -latitude,
        T ~ latitude
    ))


# # Metadata finalisation
#
# ## Split dataset
#
# Until now we have combined the per-accession metadata with the run-level
# metadata for simplicity. However, we want the two separate, so split them
# into two tables.

str(all_meta)

all_meta = all_meta %>%
    mutate(dataset=substr(oa_id, 3, 4))

all_sra = all_meta %>%
    filter(is.na(internal_data) | !internal_data) %>%
    select(sra_run, biosample, bioproject, instrument_model, library_layout, oa_id) %>%
    unique()

all_acc = all_meta %>%
    select(oa_id, dataset, ecotype_id, sample_name, biosample, species,
           latitude, longitude, elevation, locality, country, collector,
           collection_date, collection_notes, cs_number, ploidy, internal_data,
           pooled_plants) %>%
    unique()


# Did we miss any columns?
setdiff(union(colnames(all_acc), colnames(all_sra)), colnames(all_meta))
setdiff(colnames(all_meta), union(colnames(all_acc), colnames(all_sra)))

str(all_acc)
str(all_sra)

table(all_acc$dataset)


# ## SRA metadata validation
#
# This is pretty cursary, chekcing that the run IDs look legit, and there isn't
# any duplication.

sra_val = validator(
    run = grepl("(ERR|SRR|DRR)\\d+", sra_run, perl=T),
    run_unq = multiplicity(sra_run) == 1,
    bioproj = !is.na(bioproject),
    biosample = !is.na(biosample),
    liblay = library_layout %in% c("PAIRED", "SINGLE")
)

all_sra_val = confront(all_sra, sra_val)
summary(all_sra_val)
errors(all_sra_val)

all_sra_bad = all_sra %>%
    mutate(failed = get_failing_tests(all_sra_val)) %>%
    filter(failed != "")
write_tsv(all_sra_bad, "all_sra_bad.tsv")


# ## Accession metadata

acc_val = validator(
    oa_unq = multiplicity(oa_id) == 1,
    eco_unq = multiplicity(ecotype_id) == 1,
    # Some names are duplicated, but have unique ecotype IDs.
    name_unq  = (!is.na(ecotype_id)) | multiplicity(sample_name) == 1,
    latlong_ok = !(is.na(latitude) | is.na(longitude)),
    dist2land = match_country(longitude, latitude, 'dist') < 5,
    ctrymatch = match_country(longitude, latitude, 'ISO3') == country_code,
    spp = !is.na(species) | !grepl("^Arabidopsis", species)
)

all_acc_val = confront(all_acc, acc_val)
summary(all_acc_val)
errors(all_acc_val)

all_acc_bad = all_acc %>%
    mutate(failed = get_failing_tests(all_acc_val)) %>%
    filter(failed != "")

write_tsv(all_acc_bad, "all_acc_bad.tsv")

# ## Acanthophis sample2runlib file

rl2s = all_meta %>%
    transmute(
        library=ifelse(is.na(sra_run), internal_sample_name, sra_run),
        run=ifelse(is.na(bioproject), dataset, bioproject),
        sample=oa_id,
        include="Y",
        is_sra=ifelse(is.na(internal_data)| !internal_data, "Y", ""),
    )
write_tsv(rl2s, "omniath_rl2s.tsv")

# ### Sample sets

all_acc %>%
    filter(species=="Arabidopsis thaliana") %>%
    pull(oa_id) %>%
    writeLines("samplesets/Athaliana.txt")

all_acc %>%
    filter(species!="Arabidopsis thaliana") %>%
    pull(oa_id) %>%
    writeLines("samplesets/nonAthaliana.txt")

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

# # SRA size and coverage estimation

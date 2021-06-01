# # Omni-ath metadata import, validation, and merge
#
# 

library(tidyverse)
library(readxl)
library(janitor)
library(validate)
library(sp)
#install.packages(c("tidyverse", "sp", "sf", "ggplot2", "ggmap", "validate", "readxl", "janitor"))


# ## The OG 1001 Genomes
#
# NB: this data extracted from Joffery's 1kg website.

kg_sra = read_csv("source/1kg/1kg-sra-run-table.csv")
kg_acc = read_csv("source/1kg/1kg-accessions.csv")

# ## African lines
#
# From Angela Hancock's lab

af_sra = read_tsv("source/africans/african_sra-run-table-PRJEB19780.tsv")
af_acc = read_csv("source/africans/tabula-pnas.1616736114.sapp.csv", na='-')


# ## Argentinian lines
#
# From [Luciana Kasulin's paper in MolEcol](https://dx.doi.org/10.1111/mec.14107)

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


# ## Chinese lines
#
# From Yangtze river basin, see [Zou et al.
# 2017](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1378-9).


chi_sra = read_tsv("source/chinese/SraRunTable_PRJNA293798.txt")
chi_acc = read_tsv("source/chinese/supp-table-3.tsv")


# ## Catalonia
#
# https://www.pnas.org/content/115/52/E12443

cat_sra = read_tsv("source/catalonian/catalonianSraRunTable.txt")

cat_acc_orig = read_xlsx("source/catalonian/pnas.1816964115.sd01.xlsx") %>%
    clean_names()

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
    select(deme_of_origin, town, gps_coordinates, location) %>%
    filter(!is.na(town)) %>%
    mutate(gps = purrr::map(gps_coordinates, function(x) {
            g = str_split(x, ",") %>%
                unlist() %>%
                purrr::map_chr(str_trim) %>%
                as.numeric()
            data.frame(latitude=g[1], longitude=g[2])
    })) %>%
    unnest(gps) %>%
    select(-gps_coordinates)

cat_acc = cat_acc_sample %>%
    left_join(cat_acc_deme, by="deme_of_origin")


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
tib_acc = tribble(
    ~sample, ~latitude, ~longitude, ~elevation,
    "Tibet_0", 29.6903, 90.9338, 4200,
)


# # Merge

# ## Accession metadata
#

kg_acc_pre = kg_acc %>%
    select(ecotype_id, sample_name=name, latitude, longitude, locality, 
           country, collector, collection_date, cs_number)
kg_acc_pre
af_acc_pre = af_acc %>%
    mutate(collection_date = as.character(collection_year)) %>%
    select(sample_name=name, latitude=Latitude, longitude=Longitude,
           locality=Region, country=Country, collection_date)
af_acc_pre

arg_acc_pre = arg_acc %>%
    mutate(country="Argentina", locality="Patagonia") %>%
    select(sample_name=name, latitude, longitude, elevation, collector,
           country, locality)
arg_acc_pre

chi_acc_pre = chi_acc %>%
    select(sample_name=SampleID, latitude=Latitude, longitude=Longitude,
           locality=Region)
chi_acc_pre

cat_acc_pre = cat_acc %>%
    select(sample_name, latitude, longitude, locality=town)
cat_acc_pre
View(cat_acc)

tib_acc_pre = tib_acc %>%
    select(sample_name=sample, latitude, longitude, elevation)
tib_acc_pre

all_acc = bind_rows(kg_acc_pre, af_acc_pre, arg_acc_pre, chi_acc_pre,
                    cat_acc_pre, tib_acc_pre)


# ### Consistent unique IDs for all accessions
#
# Not all plants have ecotype IDs, and many names have charachers we dont' want
# in file names (non-ascii, or weird punctuation). For our general sanity, we
# want a simple ID like an ecotype ID for all samples. For samples with ecotype
# IDs, we'll just use the ecotype ID. For others, we'll make new ones.


all_acc_km = all_acc %>%
    mutate(km_id = ifelse(is.na(ecotype_id),
                          sprintf("km%04d", 1:n()),
                          sprintf("at%d", ecotype_id)))



# ## Run/SRA level metadata
#
# Here we pluck a core set of fields from each SRA table, then merge it into
# one large table for everything.


kg_sra_pre = kg_sra %>%
    select(sra_run=Run, bioproject=BioProject, ecotype_id=Ecotype,
           library_layout=LibraryLayout, instrument=Instrument) %>%
    left_join(all_acc_km %>%
        filter(!is.na(ecotype_id)) %>%
        select(ecotype_id, km_id),
        by="ecotype_id")
stopifnot(all(!is.na(kg_sra_pre$km_id)))
kg_sra_pre


# Just the name:km_id mapping for non-1kg lines
acc_just_ids = all_acc_km %>%
        filter(is.na(ecotype_id)) %>%
        select(sample_name, km_id)


af_sra_pre = af_sra %>%
    select(sra_run=run_accession, bioproject=study_accession,
           sample_name=sample_alias, library_layout,
           instrument=instrument_model) %>%
    left_join(acc_just_ids, by="sample_name")
stopifnot(all(!is.na(af_sra_pre$km_id)))
af_sra_pre


arg_sra_pre = arg_sra %>%
    select(sra_run=run_accession, bioproject=study_accession,
           sample_name=sample_alias, library_layout,
           instrument=instrument_model) %>%
    left_join(acc_just_ids, by="sample_name")
stopifnot(all(!is.na(arg_sra_pre$km_id)))
arg_sra_pre


chi_sra_pre = chi_sra %>%
    select(sra_run=run_accession, bioproject=study_accession,
           sample_name=sample_alias, library_layout,
           instrument=instrument_model) %>%
    left_join(acc_just_ids, by="sample_name")
stopifnot(all(!is.na(chi_sra_pre$km_id)))
chi_sra_pre


cat_sra_pre = cat_sra %>%
    select(sra_run=Run, bioproject=BioProject, sample_name=Sample_Name,
           library_layout=LibraryLayout, instrument=Instrument) %>%
    left_join(acc_just_ids, by="sample_name")
stopifnot(all(!is.na(cat_sra_pre$km_id)))
cat_sra_pre

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


# Merge both dataset after above "fixing". 

cat_merged = full_join(cat_sra_pre_name, cat_acc, by=c("match_name"="sample_name")) %>%
    arrange(is.na(deme), is.na(deme_of_origin), sample_name)
cat_merged

# Most samples

sra_deme = table(cat_merged$deme)
acc_deme = table(cat_merged$deme_of_origin)
length(sra_deme)
length(acc_deme)


tib_sra_pre = tib_sra %>%
    select(sra_run=Run, bioproject=BioProject, sample_name=`Sample Name`,
           library_layout=LibraryLayout, instrument=Instrument) %>%
    left_join(acc_just_ids, by="sample_name")
tib_sra_pre


all_sra = bind_rows(kg_sra_pre, af_sra_pre, arg_sra_pre, chi_sra_pre,
                    cat_sra_pre, tib_sra_pre)


View(cat_sra)

View(all_sra)
# # Validate metadata
#
# ## SRA

str(all_sra)

val = validator(
    run = !is.na(sra_run) & grepl("(ERR|SRR)\\d+", sra_run, perl=T),
    bioproj = !is.na(bioproject),
    liblay = library_layout %in% c("PAIRED", "SINGLE"),
    name_ok = !(is.na(ecotype_id) & is.na(sample_name))
)

sra_val = confront(all_sra, val)
summary(sra_val)

sra_bad = violating(all_sra, val)
sra_bad


# ## Accession metadata

multiplicity = function(x) {
    x = as.character(x)
    table(x)[x]
}

val = validator(
    name_ok = !(is.na(ecotype_id) & is.na(sample_name)),
    eco_unq = multiplicity(ecotype_id) == 1,
    # Some names are duplicated, but have unique ecotype IDs.
    name_unq  = (!is.na(ecotype_id)) | multiplicity(sample_name) == 1,
    latlong_ok = !(is.na(latitude) | is.na(longitude))
)

acc_val = confront(all_acc, val)
summary(acc_val)

acc_bad = violating(all_acc, acc_val)
acc_bad




#' # Find non-ACGTN bases in all references
#'
#' We have an issue with bcftools norm due to non-ACGTN bases appearing in
#' references. This code should determine which references this affects, so I
#' can either correct them, or update the REF column of the VCF files.


import fastaparser
import glob
import yaml
import pandas as pd

with open("../../config.yml") as fh:
    snkconf = yaml.load(fh, Loader=yaml.FullLoader)

ACGTN = set(list("ACGTN"))
nonacgtn = []
for ref, dat in snkconf["refs"].items():
    with open("../../"+dat["fasta"]) as fh:
        fas = fastaparser.Reader(fh, parse_method="quick")
        for ctg in fas:
            ctgname = ctg.header.lstrip(">").strip().split(" ")[0]
            for i, base in enumerate(ctg.sequence):
                if base.upper() not in ACGTN:
                    nonacgtn.append({
                        "ref": ref, "contig": ctgname, "pos": i+1,
                        "base": base.upper()
                    })



#' Now we have essentialy a BED file of non-ACGTN bases. Let's
#' use pandas to summarise the per-base observaitions by ref an
#' contig.

ndf = pd.DataFrame(nonacgtn)
ndf.groupby(["ref", "contig"])["contig"].count()

#' So the issues are all with TAIR10. I'll have to run the
#' bcftools plugin to update the reference.


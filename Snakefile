configfile: "config.yml"
import acanthophis
acanthophis.populate_metadata(config)


include: acanthophis.rules.base
include: "rules/sra.rules"
include: "rules/reads.rules"
include: "rules/qualimap.rules"
include: acanthophis.rules.align
include: acanthophis.rules.varcall
include: acanthophis.rules.multiqc
include: acanthophis.rules.kraken
#include: acanthophis.rules.variantannotation
include: "rules/taxonid.rules"


# for running on taco etc
localrules: bwaidx, mergebam_samp, bamidx, bamstat_sample, qualimap_set,  bcfnorm, bcffilter, multiqc_bamstats, multiqc_kraken, split_pairs_r12, split_pairs_il, qcreads_se, read_count_librun_indiv, samplefastqfile, multiqc_rawreads, multiqc_samplereads
#localrules: bwaidx, mergebam_samp, bam_markdups_sort, bamidx, bamstat_sample, qualimap_set,  bcfnorm, bcffilter, bcfmerge, multiqc_bamstats, multiqc_kraken, split_pairs_r12, split_pairs_il, qc_reads_se, read_count_librun_raw, samplefastqfile, multiqc_rawreads, multiqc_samplereads

rule all:
    input:
        rules.reads.input,
        rules.align.input,
        rules.varcall.input,
        rules.multiqc.input,
        rules.all_kraken.input,
#        rules.all_snpeff.input,


align_rule_samples = set()
for sset in config["align"]["samplesets"]:
    for s in config['SAMPLESETS'][sset]:
        align_rule_samples.add(s)

localrules: align_samples_idx
rule align_samples_idx:
    input:
        expand("data/alignments/samples/{aligner}/{ref}/{sample}.bam",
               ref=config["align"]["refs"],
               aligner=config["align"]["aligners"],
               sample=align_rule_samples),
        expand("data/alignments/samples/{aligner}/{ref}/{sample}.bam.bai",
               ref=config["align"]["refs"],
               aligner=config["align"]["aligners"],
               sample=align_rule_samples),
        expand("data/alignments/bamstats/sample/{aligner}~{ref}~{sample}.samtools.stats",
               ref=config["align"]["refs"],
               aligner=config["align"]["aligners"],
               sample=align_rule_samples),

localrules: align_persampleset
rule align_persampleset:
    input:
        [ expand(path, ref=c["refs"], aligner=c["aligners"], sample=config["SAMPLESETS"][s])
          for s, c in config["align"]["persampleset"].items()
          for path in (
              "data/alignments/samples/{aligner}/{ref}/{sample}.bam",
              "data/alignments/samples/{aligner}/{ref}/{sample}.bam.bai",
              "data/alignments/bamstats/sample/{aligner}~{ref}~{sample}.samtools.stats",
          )]

rule difflines_qualimap:
    input:
        expand("data/alignments/qualimap/sampleset/{aligner}~{ref}~{sampleset}/", 
               ref=config["align"]["refs"],
               aligner=config["align"]["aligners"],
               sampleset="difflines")


rule fp_taxonid:
    input:
        [ expand("data/taxonid/{sampleset}/kraken/{krakendb}/{sample}_unclassified.fastq.gz",
                 sampleset=sampleset, krakendb=dat["kraken"]["dbs"], sample=config["SAMPLESETS"][sampleset])
          for sampleset, dat in config["taxonid"]["samplesets"].items()],
        [ expand("data/taxonid/{sampleset}/diamond/{krakendb}/{diamonddb}~{sample}.daa",
                 sampleset=sampleset, krakendb=dat["kraken"]["dbs"], diamonddb=dat["diamond"]["dbs"], sample=config["SAMPLESETS"][sampleset])
          for sampleset, dat in config["taxonid"]["samplesets"].items()],

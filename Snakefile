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

rule difflines_qualimap:
    input:
        expand("data/alignments/qualimap/sampleset/{aligner}~{ref}~{sampleset}/", 
               ref=config["align"]["refs"],
               aligner=config["align"]["aligners"],
               sampleset="difflines")



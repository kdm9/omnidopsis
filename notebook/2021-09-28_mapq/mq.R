library(tidyverse)
library(foreach)
library(doParallel)
registerDoParallel(16)

refs = c("at6137", "at6923", "at6929", "at7143", "at8285", "at9104", "at9336",
         "at9503", "at9578", "at9744", "at9762", "at9806", "at9830", "at9847",
         "at9852", "at9879", "at9883", "at9900", "cdm0", "jiao_an1",
         "jiao_c24", "jiao_cvi", "jiao_eri", "jiao_kyo", "jiao_ler",
         "jiao_sha", "TAIR10")
samples = readLines("rawdata/metadata/samplesets/difflines.txt")

mq = foreach(ref=refs, .combine=rbind)%:%
foreach(sample=samples, .combine=rbind) %dopar% {
    bampath = sprintf("data/alignments/samples/bwa/%s/%s.bam", ref, sample)
    ctgmq = system(sprintf("samtools view -F 1796 %s | cut -f 3,5", bampath),
                   intern=T) %>%
        read.table(text=., header=F, sep="\t", col.names=c("contig", "mq")) %>%
        group_by(contig) %>%
        summarise(mq=mean(mq), nreads=n())
    cbind(ref=ref, sample=sample, ctgmq)
}

write_tsv(mq, "mq.tsv")

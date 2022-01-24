set -x
mkdir -p data/filtered
for ref in 9744 9852
do
    bcftools filter \
        --IndelGap 10 --SnpGap 10 \
        -i 'F_missing < 0.2 & MAC >= 3 & QUAL > 100 & INFO/DP > 200 & INFO/DP < 20000'  \
        -Oz -o data/filtered/freebayes~bwa~at${ref}~filter-miss20-mac3-qual100-dp200-20000.vcf.gz \
        data/raw/freebayes~bwa~at${ref}~difflines~filtered-default.vcf.gz
done


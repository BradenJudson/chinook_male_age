#!/bin/bash

cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/chinook_RADseq
source ~/miniconda3/etc/profile.d/conda.sh
conda activate ../hdgbs_env

# Original VCF (genome wide).
FULLVCF="./08c_filtering_paralogs/snps_maf001_singletons_sub.recode.vcf.gz"

CHRS=`cat <(bcftools query -f "%CHROM\n" ./08c_filtering_paralogs/snps_maf001_singletons_sub.recode.vcf.gz) | sort -u`

conda activate ./plink_env

for chromosome in $CHRS
do
        plink --vcf "$FULLVCF" --chr "$chromosome" --recode --out ./10_popLDs/"$chromosome"_LDs --aec --r2 --inter-chr --ld-window-r2 0
done


#!/bin/bash

cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/chinook_RADseq
source ~/miniconda3/etc/profile.d/conda.sh
conda activate ./plink_env

# Original VCF (genome wide).
FULLVCF="./08c_filtering_paralogs/snps_maf001_singletons_sub.recode.vcf.gz"

# Index file for filtering. Only required once.
#tabix "$FULLVCF"

OTS17="NC_056445.1"
OTS18="NC_056446.1"
OTS30="NC_056458.1"

# Calculate LD for chromosomes of interest.
plink --vcf "$FULLVCF" --chr "$OTS17" --recode --out ./09_Ots17/Ots17 --aec --r2 --inter-chr --ld-window-r2 0
plink --vcf "$FULLVCF" --chr "$OTS18" --recode --out ./09b_Ots18/Ots18 --aec --r2 --inter-chr --ld-window-r2 0
plink --vcf "$FULLVCF" --chr "$OTS30" --recode --out ./09c_Ots30/Ots30 --aec --r2 --inter-chr --ld-window-r2 0


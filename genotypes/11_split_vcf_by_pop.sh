#!/bin/bash

cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/chinook_RADseq
source ~/miniconda3/etc/profile.d/conda.sh
conda activate ../hdgbs_env

FULLVCF="./08c_filtering_paralogs/snps_maf001_singletons_sub.recode.vcf"
INFO="./01_info_files/sample_pops_n707.txt"
POPLIST=`awk '{print $1}' ./01_info_files/sample_pops_n707.txt | sort | uniq`

for POP in $POPLIST
do
        echo $POP

	# Extract sample IDs that correspond to each population and save them as a text file.
        awk -v population="$POP" '$1 == population {print $2}' "$INFO" > ./01_info_files/"$POP"_indvs.txt

	# Subset the full VCF to only individuals within the population of interest.
	# Exclude monomorphic alleles (i.e., minor allele count = 0). 
        bcftools view --samples-file ./01_info_files/"$POP"_indvs.txt "$FULLVCF" \
                | bcftools view -e "MAC == 0" --output-type v --output ./08d_pop_vcfs/"$POP".vcf

	# Take the above pop-specific VCF and conduct linkage pruning.
        bcftools +prune --max 0.2 -w 10000 ./08d_pop_vcfs/"$POP".vcf  --output-type v --output ./08d_pop_vcfs/"$POP"_LDprune.vcf
done
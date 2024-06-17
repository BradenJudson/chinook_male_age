#!/bin/bash

cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/chinook_RADseq
source ~/miniconda3/etc/profile.d/conda.sh
conda activate ./plink_env

POPLIST=`awk '{print $1}' ./01_info_files/sample_pops_n707.txt | sort | uniq`
fastPHASE="../software/fastphase/fastPHASE"

for POP in $POPLIST
do
        conda activate ./plink_env

        # These need to be bgzipped and index with tabix beforehand.
        VCF=./08d_pop_vcfs/"$POP".vcf.gz
	
	# Isolate Ots17 and get pairwise R2 between all SNPs.
        plink --vcf "$VCF" --chr "NC_056445.1" --recode --out ./11_popOts17/"$POP"_Ots17 --aec --r2 --inter-chr --ld-window-r2 0
	
	# Extract only highly linked (R2 > 0.3) SNPs and extract their positions along the chromosome.
        awk '{ if($7 > 0.30) { print $2,$5}}' ./11_popOts17/"$POP"_Ots17.ld > ./11_popOts17/"$POP"_SNPs.txt

	# Make a list of all unique SNPs in the highly-linked set isolated above and add a column with the chromosome name for further filtering.
        sort -u <(cat ./11_popOts17/"$POP"_SNPs.txt | grep -v "^BP"  | cut -d ' ' -f1  ./11_popOts17/"$POP"_SNPs.txt | grep -v "^BP"  | cut -d ' ' -f2) | \
               awk '{gsub(/ /,",",$0); print "NC_056445.1\t"$0}' > ./11_popOts17/"$POP"_highLD_SNPpos.txt	
	
	# For each population, extract only the highly linked SNPs identified above.
        bcftools view --regions-file ./11_popOts17/"$POP"_highLD_SNPpos.txt --output ./11_popOts17/"$POP"_HLD.vcf --output-type v "$VCF"

	# Convert the highly linked SNPs for each population to the fastPHASE format using plink.
        plink --vcf ./11_popOts17/"$POP"_HLD.vcf --recode fastphase --out ./11_popOts17/"$POP"_fastphase --aec

	# Phase the corresponding genotypes and write.
        ./"$fastPHASE" ./11_popOts17/"$POP"_fastphase.chr-NC_056445.1.recode.phase.inp -o"$POP"_phased

done
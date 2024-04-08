#!/bin/bash


#SBATCH --job-name=radseq_cutadapt
#SBATCH --partition=standard
#SBATCH --mem-per-cpu=5G
#SBATCH --mail-user=Braden.Judson@dfo-mpo.gc.ca
#SBATCH --time=16:00:00
#SBATCH --cpus-per-task=64
#SBATCH --account=grdi_genarcc


source ~/.profile


NCPU=64

# Load GNU parallel.
. ssmuse-sh -x main/opt/parallel/parallel-202109022

# Export for parallelization.
export cutadapt

#/gpfs/fs7/grdi/genarcc/wp3/judsonb/radseq_env/bin/cutadapt

ls -1 02a_raw_reads/*_R1.fastq.gz |
         sed 's/_R1\.fastq\.gz//g' |
         sed 's/02a\_raw_reads\///g' |
         parallel -j $NCPU cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
         -g AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
         -o 02b_reads_trimmed/{}"_R1.fastq.gz" \
         -p 02b_reads_trimmed/{}"_R2.fastq.gz" \
         -e 0.2 \
         -m 50 \
            02a_raw_reads/{}_R1.fastq.gz \
            02a_raw_reads/{}_R2.fastq.gz

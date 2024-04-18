#!/bin/bash

#SBATCH --job-name=radseq_gstacks
#SBATCH --partition=standard
#SBATCH --mail-user=Braden.Judson@dfo-mpo.gc.ca
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=5G
#SBATCH --account=grdi_genarcc

# Source conda environment.
source ~/.profile

STACKS_FOLDER="06_stacks"
SAMPLE_FOLDER="05_alignment"
INFO_FILES="01_info_files"
POP_MAP="population_map.txt"

gstacks -I "$SAMPLE_FOLDER" -S ".1.sorted.bam" \
        -M "$INFO_FILES"/"$POP_MAP" \
        -O "$STACKS_FOLDER" \
        --max-clipped 0.1
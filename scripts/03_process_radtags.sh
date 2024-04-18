#!/bin/bash


#SBATCH --job-name=radseq_pradtags
#SBATCH --partition=standard
#SBATCH --mem-per-cpu=5G
#SBATCH --ntasks=1
#SBATCH --mail-user=Braden.Judson@dfo-mpo.gc.ca
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=64
#SBATCH --account=grdi_genarcc


source ~/.profile

ENZYME1="PstI"
ENZYME2="MspI"
TRIM_LENGTH=130
NCPU=64


# Load GNU parallel.
. ssmuse-sh -x main/opt/parallel/parallel-202109022

export stacks
export process_radtags

INFO_FILES="01_info_files"

# Pipe works such that the lane info strings are read as the third argument in the utility script.

cat $INFO_FILES/lane_info.txt |
        parallel -j $NCPU ./00_scripts/utility_scripts/process_radtags_2_enzyme_pe.sh "$TRIM_LENGTH" "$ENZYME1" "$ENZYME2"

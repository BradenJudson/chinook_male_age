#!/bin/bash


ENZYME1="PstI"
ENZYME2="MspI"
TRIM_LENGTH=130
NCPU=64

source ~/.bashrc

# Load GNU parallel.
. ssmuse-sh -x main/opt/parallel/parallel-202109022

export stacks
export process_radtags

INFO_FILES="01_info_files"

cat $INFO_FILES/lane_info.txt |
        parallel -j $NCPU ./00_scripts/utility_scripts/process_radtags_1_enzyme_pe.sh "$TRIM_LENGTH" "$ENZYME1" "$ENZYME2"




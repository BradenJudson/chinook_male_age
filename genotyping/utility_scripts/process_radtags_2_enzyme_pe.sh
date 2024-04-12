#!/bin/bash


# Source conda environment.
source ~/miniconda3/etc/profile.d/conda.sh
cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/chinook_RADseq
conda activate radseq_env


TRIM_LENGTH=$1
ENZYME1=$2
ENZYME2=$3
LANE=$4
INFO_FILES="01_info_files"
TEMPBARCODES=.temp.${LANE}.barcodes


grep -vE "^$" "$INFO_FILES"/chRADseq_sampleinfo_stacks.csv | \
    grep "$LANE" | \
    grep -v "Barcode" | \
    cut -d ',' -f 2 > "$INFO_FILES"/"$TEMPBARCODES"

process_radtags -i gzfastq \
        -1 02c_reads_trimmed_concat/"$LANE"".fastq.gz" \
        -2 02c_reads_trimmed_concat/$(echo "$LANE" | perl -pe 's/_R1/_R2/')".fastq.gz" \
        -o 03_samples/"$LANE" \
        -b "$INFO_FILES"/"$TEMPBARCODES" \
        -c -r -t "$TRIM_LENGTH" \
        -q \
        -s 0 \
        -E phred33 \
        --barcode-dist-1 2 \
        --renz_1 "$ENZYME1" \
        --renz_2 "$ENZYME2"
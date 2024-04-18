#!/bin/bash

cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/chinook_RADseq/02c_reads_trimmed_concat

for file in *.gz
do
        sample="$file"
        nreads=$(zcat "$file" | echo $((`wc -l`/4)));
        echo ${sample} ${nreads}
done > ../01_info_files/sample_reads.txt


#!/bin/bash


#SBATCH --job-name=merge_plates
#SBATCH --partition=standard
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=8
#SBATCH --mail-user=Braden.Judson@dfo-mpo.gc.ca
#SBATCH --time=01:00:00
#SBATCH --account=grdi_genarcc


cd /gpfs/fs7/grdi/genarcc/wp3/judsonb/chinook_RADseq/02_reads


# Plate 1, forward:
cat NS.LH00487_0020.005.D705---B502.Thealy-BC1-96-P1_R1.fastq.gz \
    NS.LH00487_0020.006.D705---B502.Thealy-BC1-96-P1_R1.fastq.gz \
    NS.LH00487_0020.007.D705---B502.Thealy-BC1-96-P1_R1.fastq.gz > ../02c_reads_trimmed_concat/trimmed_radseq_P1_R1.fastq.gz

# Plate 1, reverse:
cat NS.LH00487_0020.005.D705---B502.Thealy-BC1-96-P1_R2.fastq.gz \
    NS.LH00487_0020.006.D705---B502.Thealy-BC1-96-P1_R2.fastq.gz \
    NS.LH00487_0020.007.D705---B502.Thealy-BC1-96-P1_R2.fastq.gz > ../02c_reads_trimmed_concat/trimmed_radseq_P1_R2.fastq.gz

# Plate 2, forward:
cat NS.LH00487_0020.005.D705---B503.Thealy-BC97-192-P2_R1.fastq.gz \
    NS.LH00487_0020.006.D705---B503.Thealy-BC97-192-P2_R1.fastq.gz \
    NS.LH00487_0020.007.D705---B503.Thealy-BC97-192-P2_R1.fastq.gz > ../02c_reads_trimmed_concat/trimmed_radseq_P2_R1.fastq.gz

# Plate 2, reverse:
cat NS.LH00487_0020.005.D705---B503.Thealy-BC97-192-P2_R2.fastq.gz \
    NS.LH00487_0020.006.D705---B503.Thealy-BC97-192-P2_R2.fastq.gz \
    NS.LH00487_0020.007.D705---B503.Thealy-BC97-192-P2_R2.fastq.gz > ../02c_reads_trimmed_concat/trimmed_radseq_P2_R2.fastq.gz

# Plate 3, forward:
cat NS.LH00487_0020.005.D705---B504.Thealy-BC193-288-P3_R1.fastq.gz \
    NS.LH00487_0020.006.D705---B504.Thealy-BC193-288-P3_R1.fastq.gz \
    NS.LH00487_0020.007.D705---B504.Thealy-BC193-288-P3_R1.fastq.gz > ../02c_reads_trimmed_concat/trimmed_radseq_P3_R1.fastq.gz

# Plate 3, reverse:
cat NS.LH00487_0020.005.D705---B504.Thealy-BC193-288-P3_R2.fastq.gz \
    NS.LH00487_0020.006.D705---B504.Thealy-BC193-288-P3_R2.fastq.gz \
    NS.LH00487_0020.007.D705---B504.Thealy-BC193-288-P3_R2.fastq.gz > ../02c_reads_trimmed_concat/trimmed_radseq_P3_R2.fastq.gz

# Plate 4, forward:
cat NS.LH00487_0020.005.D706---B502.Thealy-BC289-384-P4_R1.fastq.gz \
    NS.LH00487_0020.006.D706---B502.Thealy-BC289-384-P4_R1.fastq.gz \
    NS.LH00487_0020.007.D706---B502.Thealy-BC289-384-P4_R1.fastq.gz > ../02c_reads_trimmed_concat/trimmed_radseq_P4_R1.fastq.gz

# Plate 4, reverse:
cat NS.LH00487_0020.005.D706---B502.Thealy-BC289-384-P4_R2.fastq.gz \
    NS.LH00487_0020.006.D706---B502.Thealy-BC289-384-P4_R2.fastq.gz \
    NS.LH00487_0020.007.D706---B502.Thealy-BC289-384-P4_R2.fastq.gz > ../02c_reads_trimmed_concat/trimmed_radseq_P4_R2.fastq.gz


# Plate 5, forward:
cat NS.LH00487_0020.007.D706---B503.Thealy-BC1-96-P5_R1.fastq.gz \
    NS.LH00487_0020.006.D706---B503.Thealy-BC1-96-P5_R1.fastq.gz \
    NS.LH00487_0020.005.D706---B503.Thealy-BC1-96-P5_R1.fastq.gz > ../02c_reads_trimmed_concat/trimmed_radseq_P5_R1.fastq.gz

# Plate 5, reverse:
cat NS.LH00487_0020.007.D706---B503.Thealy-BC1-96-P5_R2.fastq.gz \
    NS.LH00487_0020.006.D706---B503.Thealy-BC1-96-P5_R2.fastq.gz \
    NS.LH00487_0020.005.D706---B503.Thealy-BC1-96-P5_R2.fastq.gz > ../02c_reads_trimmed_concat/trimmed_radseq_P5_R2.fastq.gz

# Plate 6, forward:
cat NS.LH00487_0020.005.D706---B504.Thealy-BC97-192-P6_R1.fastq.gz \
    NS.LH00487_0020.006.D706---B504.Thealy-BC97-192-P6_R1.fastq.gz \
    NS.LH00487_0020.007.D706---B504.Thealy-BC97-192-P6_R1.fastq.gz > ../02c_reads_trimmed_concat/trimmed_radseq_P6_R1.fastq.gz

# Plate 6, reverse:
cat NS.LH00487_0020.005.D706---B504.Thealy-BC97-192-P6_R2.fastq.gz \
    NS.LH00487_0020.006.D706---B504.Thealy-BC97-192-P6_R2.fastq.gz \
    NS.LH00487_0020.007.D706---B504.Thealy-BC97-192-P6_R2.fastq.gz > ../02c_reads_trimmed_concat/trimmed_radseq_P6_R2.fastq.gz

# Plate 7, forward:
cat NS.LH00487_0020.005.D707---B503.Thealy-BC193-288-P7_R1.fastq.gz \
    NS.LH00487_0020.006.D707---B503.Thealy-BC193-288-P7_R1.fastq.gz \
    NS.LH00487_0020.007.D707---B503.Thealy-BC193-288-P7_R1.fastq.gz > ../02c_reads_trimmed_concat/trimmed_radseq_P7_R1.fastq.gz

# Plate 7, reverse:
cat NS.LH00487_0020.005.D707---B503.Thealy-BC193-288-P7_R2.fastq.gz \
    NS.LH00487_0020.006.D707---B503.Thealy-BC193-288-P7_R2.fastq.gz \
    NS.LH00487_0020.007.D707---B503.Thealy-BC193-288-P7_R2.fastq.gz > ../02c_reads_trimmed_concat/trimmed_radseq_P7_R2.fastq.gz

# Plate 8, forward:
cat NS.LH00487_0020.005.D707---B504.Thealy-BC289-384-P8_R1.fastq.gz \
    NS.LH00487_0020.006.D707---B504.Thealy-BC289-384-P8_R1.fastq.gz \
    NS.LH00487_0020.007.D707---B504.Thealy-BC289-384-P8_R1.fastq.gz > ../02c_reads_trimmed_concat/trimmed_radseq_P8_R1.fastq.gz

# Plate 8, reverse:
cat NS.LH00487_0020.005.D707---B504.Thealy-BC289-384-P8_R2.fastq.gz \
    NS.LH00487_0020.006.D707---B504.Thealy-BC289-384-P8_R2.fastq.gz \
    NS.LH00487_0020.007.D707---B504.Thealy-BC289-384-P8_R2.fastq.gz > ../02c_reads_trimmed_concat/trimmed_radseq_P8_R2.fastq.gz





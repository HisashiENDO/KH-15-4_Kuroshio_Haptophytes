#!/bin/sh

# This is the main scropt to sequence trimming and quality filtering for QIIME2 pipeline
# All the fastq files used were merged as "merge.fastq"

# set directory

module load seqkit/0.14.0
module load cutadapt/2.5
module load fastx_toolkit/0.0.14


#### Quality filtering ####
# Cut primer sequences (remove reverse first, then forward and do size selection)
## A1 and trP1 adapters were already trimmed off
## Reverse: CCAGGGATGTTTTCACTGATC (GATCAGTGAAAACATCCCTGG)
cutadapt -e 0.3 -a CCAGGGATGTTTTCACTGATC --discard-untrimmed merge.fastq  > merge.fastq.remr.fastq
## Forward: AGCCGCGGTAATTCCA
cutadapt -e 0.3 -g AGCCGCGGTAATTCCA -m 300 -M 500 --discard-untrimmed merge.fastq.remr.fastq > merge.fastq.remrf.fastq
# -g forward, -a reverse, -m min length, -M max length, --discard-untrimmed Discard reads that do not contain an adapter.

# Trim off the region after 250 bp to keep high quality reads
fastx_trimmer -l 250 -i merge.fastq.remrf.fastq -o merge.fastq.remrf.trim250.fastq

# Trim low quality reads
fastq_quality_filter -v -i merge.fastq.remrf.trim250.fastq -o merge.fastq.remrf.trim250.filt.fastq -q 25 -p 90 -Q 33

# Stats
fastx_quality_stats -i merge.fastq.remrf.trim250.filt.fastq -o merge.fastq.remrf.trim250.filt.fastq.text -Q33
fastq_quality_boxplot_graph.sh -i merge.fastq.remrf.trim250.filt.fastq.text -o merge.fastq.remrf.trim250.filt.fastq.png -t "merge.fastq.qualfiltered"

# Show statistics for each library
seqkit grep -nirp S01_* merge.fastq.remrf.trim250.filt.fastq | seqkit stats -a

# You can change "S01_" to something

# Separate as each library
sh separate_seqs.sh merge.fastq.remrf.trim250.filt.fastq
####

module purge

echo "Quality filter finished"

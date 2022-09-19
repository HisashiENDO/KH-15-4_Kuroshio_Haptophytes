#!/bin/sh

# Separete merged fastq file into each samples by the sample ids such as S01, S02,...

FASTQ="$1"

# set directory

module load seqkit/0.14.0

mkdir curated_fastqs
out_dir="curated_fastqs/"

for i in 01 02 03 04 05 07 08 09 10 11 12 14 15 17 18 20 21 23 24 25 26 27 28 29 30 31 32 34 35 37 38 40 41 42 43 44 45 46 47 49 50 51 52 53 54 55
  do
    seqkit grep -nirp S${i}_* ${FASTQ} > ${out_dir}$i.curated.fastq
  done

module purge

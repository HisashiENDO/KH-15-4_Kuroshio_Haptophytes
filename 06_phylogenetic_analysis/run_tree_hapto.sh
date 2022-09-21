#!/bin/bash
source /etc/profile.d/modules.sh

<< COMMENTOUT
"""
# Make LM tree using curated sequences of haptophyte 18S rRNA gene
"""
COMMENTOUT

# set directory

M1=ASV280_ref63_merged.fasta

## Alignment ###
module load mafft/7.481

# Align reference polB
mafft-linsi --thread 12 ${M1} > ${M1}.msa

module purge

## Build Tree ###
module load iqtree/1.6.12 # The latest version iqtree/2.1.2 has different usage

# Select best-fit model for alignment
iqtree -s ${M1}.msa -m TESTONLY -nt	12 -pre ${M1}.msa.model.select

# -> TIM3+F+I+G4 was chosen as the best-fit model

# Tree construction
iqtree -s ${M1}.msa -nt 16 -mem 100G -m TIM3+F+I+G4 -bb 1000 -pre ${M1}.msa
# -bb 1000 Ultrafast bootstrap

module purge

# Tree construction
# qsub -q cdb -l select=1:ncpus=12:mem=100gb run_tree_hapto.sh

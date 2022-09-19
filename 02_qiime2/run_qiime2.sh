#!/bin/sh

## This script aims to make ASV profile for KH15-4 curated sequencing data
## Prepare high-quality and length unified fastq file by the pipeline in the directory "01_sequence_processing"

## The general procedures as follows
#1. Importing sequence data
#2. DADA2 denoising, defining ASVs, and making table
#3. Remove sigleton
#4. Rarefy ASVs
 # ->Then make visual sumamries for the qiime2 view
#5. Make trained classifer for the taxonomic analysis
#6. Classify ASVs into taxonomy
#7. Make bar plot of summarizing taxonomy levels
#8. Ordination based on Unifrac distance etc.
#9. Output ASV and taxonomy tables as tsv format

## QIIME2 "Moving Picture Tutrial" is available below.
# https://docs.qiime2.org/2020.8/tutorials/moving-pictures/
# Japanese private site: https://qiita.com/keisuke-ota/items/6399b2f2f7459cd9e418

# Activate QIIME2
# source activate qiime2-2020.2

# Load module
module load qiime2/2020.2
module load Python/3.7.5

# set directory

#this file can be gerenated with a script (0 MANIFEST), but manuall creation is recommended.
# tutorial: https://docs.qiime2.org/2018.11/tutorials/importing/
MANIFEST=qiime_manifest_merge.txt

#This file is needed for a barplot and other plots and has to be cre ated manually
# tutorial: https://docs.qiime2.org/2018.11/tutorials/metadata/
METADATA_FILE=sample_metadata.tsv


<< COMMENTOUT

##### Import sequences as qiime2 "qza" format #####
echo "Importing sequence data using $MANIFEST"
qiime tools import \
--type 'SampleData[SequencesWithQuality]' \
--input-path $MANIFEST \
--output-path seqs.qza \
--input-format SingleEndFastqManifestPhred33

##### Denoizng by DADA2 #####
echo "DADA2 denoising, defining ASVs, and making table"
qiime dada2 denoise-pyro \
  --p-trim-left 0 \
  --p-trunc-len 0 \
  --i-demultiplexed-seqs seqs.qza \
  --p-max-ee 2 \
  --p-chimera-method 'consensus' \
  --p-n-threads 12 \
  --o-table seqs.table.qza \
  --o-representative-sequences seqs.repsequence.qza \
  --o-denoising-stats seqs.dada2_denoising-stats.qza

# Convert to qzv format
echo "Converting DADA2 results to qzv format"
## Make visualization qza format
qiime metadata tabulate \
  --m-input-file seqs.dada2_denoising-stats.qza \
  --o-visualization seqs.dada2_denoising-stats.qzv
## Make ASV table
qiime metadata tabulate \
  --m-input-file seqs.table.qza \
  --o-visualization seqs.table.qzv

## Remove singletons
echo "Filtering rare ASVs"
qiime feature-table filter-features \
  --i-table seqs.table.qza \
  --p-min-frequency 10 \
  --p-min-samples 2 \
  --o-filtered-table seqs.table.remrare.qza
  #-p-min-frequency <- can be used to remove OTUs with less than 2 reads
  # --p-min-samples 10 <- can be used to remove OTUs that were found in less than 10 samples
qiime metadata tabulate \
  --m-input-file seqs.table.remrare.qza \
  --o-visualization seqs.table.remrare.qzv

###
COMMENTOUT


"""
# It is shown that the minimum sample has 13533 reads after the DADA2 denoising and remove rares
"""

## Rarefy table the samples that has minimum reads without replacement
echo "Rarefy ASV table"
qiime feature-table rarefy \
  --i-table seqs.table.remrare.qza \
  --p-sampling-depth 13533 \
  --p-with-replacement False \
  --o-rarefied-table seqs.table.rarefy.qza
## Make ASV table
qiime metadata tabulate \
  --m-input-file seqs.table.rarefy.qza \
  --o-visualization seqs.table.rarefy.qzv


##### Summarize feature-table #####
echo "Make visual sumamries of the data"
qiime feature-table summarize \
  --i-table seqs.table.rarefy.qza \
  --o-visualization seqs.table.rarefy.sum.qzv
  #--m-sample-metadata-file sample-metadata.tsv
# Make representative ASV sequence table with statistics
qiime feature-table tabulate-seqs \
  --i-data seqs.repsequence.qza \
  --o-visualization seqs.repsequence.qzv


##### Taxonomic analysis #####
# Full length silva dataset: https://docs.qiime2.org/2020.8/data-resources/
echo "Classify taxonomy by q2-feature-classifier plugin"
# Load taxonomic data
#wget "https://data.qiime2.org/2020.8/common/silva-138-99-seqs.qza"
#wget "https://data.qiime2.org/2020.8/common/silva-138-99-tax.qza"
##! If you use aptmp, you can link to "/aptmp/endo/Data/silva138/qiime2/"
# not included in qiime2 environent of own PC

# Define taxonomy data
Silva138_seqs=/aptmp/endo/Data/silva138/qiime2/silva-138-99-seqs.qza
Silva138_tax=/aptmp/endo/Data/silva138/qiime2/silva-138-99-tax.qza

# Extract target seqience region to the training
# https://docs.qiime2.org/2020.8/tutorials/feature-classifier/
qiime feature-classifier extract-reads \
  --i-sequences $Silva138_seqs \
  --p-f-primer AGCCGCGGTAATTCCA \
  --p-r-primer GATCAGTGAAAACATCCCTGG \
  --p-trunc-len 316 \
  --p-min-length 100 \
  --p-max-length 0 \
  --o-reads silva-138-99-seqs.curated.qza \
  --p-n-jobs 12

# Train references targeting the target sequences that was sequenced
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads silva-138-99-seqs.curated.qza \
  --i-reference-taxonomy $Silva138_tax \
  --o-classifier silva-138-99-tax.classifier.qza

# Classify taxonomy
qiime feature-classifier classify-sklearn \
  --i-classifier silva-138-99-tax.classifier.qza \
  --i-reads seqs.repsequence.qza \
  --o-classification seqs.repsequence.taxonomy.qza

# Output visualization
qiime metadata tabulate \
  --m-input-file seqs.repsequence.taxonomy.qza \
  --o-visualization seqs.repsequence.taxonomy.qzv

# Bar plot
qiime taxa barplot \
  --i-table seqs.table.rarefy.qza \
  --i-taxonomy seqs.repsequence.taxonomy.qza \
  --m-metadata-file $METADATA_FILE \
  --o-visualization seqs.taxa-bar-plots.qzv


##### Ordination based on Unifrac distance #####
# generate tree for the diversity analysis
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences seqs.repsequence.qza \
  --o-alignment seqs.repsequence.aligned.qza \
  --o-masked-alignment seqs.repsequence.masked-aligned.qza \
  --o-tree seqs.repsequence.unrooted-tree.qza \
  --p-n-threads 12 \
  --o-rooted-tree seqs.repsequence.rooted-tree.qza

# Calcurate beta-diversity
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny seqs.repsequence.rooted-tree.qza \
  --i-table seqs.table.rarefy.qza \
  --p-sampling-depth 10000 \
  --m-metadata-file $METADATA_FILE \
  --p-n-jobs 12 \
  --output-dir core-metrics-results
# --p-sampling-depth: set at the minimum sample


##### output results #####
## Make ASV table for curated table
qiime metadata tabulate \
  --m-input-file seqs.table.rarefy.qza \
  --o-visualization seqs.table.rarefy.tabulate.qzv
  #Same as "seqs.table.rarefy.qzv"

# Output data table as a tsv format
mkdir general_outputs
qiime tools extract \
  --input-path seqs.table.rarefy.tabulate.qzv \
  --output-path general_outputs

# Output taxonomy table as a tsv format
qiime tools extract \
  --input-path seqs.repsequence.taxonomy.qzv \
  --output-path general_outputs

# Output representative sequences as a tsv format
qiime tools extract \
  --input-path seqs.repsequence.qzv \
  --output-path general_outputs

echo "pipeline has finished"

###

# Execute with
# qsub -q cdb -l select=1:ncpus=12:mem=200gb run_qiime2.sh

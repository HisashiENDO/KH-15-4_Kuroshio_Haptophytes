#################### KH15-4 QIIME2 ASV table analysis ######################

#Import library
library(vegan)
library(gplots)
library(ggplot2)
library(gcookbook)
library(RColorBrewer)
library(MASS)
library(tidyverse)
library(pryr)
library(gridExtra) # Use for arrange several figures
library(otuSummary) # include "matrixConvert" to convert matrixt to df
library(ggrepel) # For add graph lavel
library(khroma) # For making color pallete
library(ggExtra) # Add density plot to the ggplot
library(Hmisc) # To calcurate pairwise r values 'as a matrix"

## Set directory
print("set working directory")

#load related data
asv_table <- read_tsv("asv_table_edit.tsv", col_names = TRUE)
copy_18S <- read_tsv("18SqPCR_for_R.table", col_names = TRUE)
env_table <- read_tsv("station_phys_chem_for_R.table", col_names = TRUE)
sample_group <- read_tsv("sample_groups_PCoA1_for_R.table", col_names = TRUE) #including sample lavel


################# Data purification #########################

# Remove a ASV that was not classified as haptophyte (checked by ETE3 tree and QIIME2 classification result)
asv_table_curate <- asv_table[, colnames(asv_table) != "e39fb07cd46053f5e6a2b93fa2f8cc12"]
# 1 sequences were removed and the output file should have col size of 438 (437 ASVs)

# Extract curated ASV IDs for the Tree analysis
curated_ASVs <- as.data.frame(colnames(asv_table_curate)[2:c(ncol(asv_table_curate))])
colnames(curated_ASVs) <- "ASVs"
write_tsv(curated_ASVs, "curated_ASVs.txt", col_names = FALSE)

##################################################################
#
#

################# Prepare key data for the diversity analyses ##################

# Make merged df
merge_df <- env_table %>% left_join(copy_18S, by= "Site", copy = FALSE) %>%
  left_join(sample_group, by= "Site", copy = FALSE) %>%
  left_join(asv_table_curate, by= "Site", copy = FALSE)

# Write data
write_tsv(merge_df, "merge_df.tsv")

#Estimate abundance and diversity
coverage <- as.data.frame(rowSums(merge_df[,21:ncol(merge_df)]))
richness <- as.data.frame(specnumber(merge_df[,21:ncol(merge_df)], MARGIN = 1))
shannon <- as.data.frame(diversity(merge_df[,21:ncol(merge_df)], MARGIN = 1, index="shannon", base = exp(1)))
evenness <- shannon/log(richness)

# Merge diversity indices with the environment
df_env_div <- merge_df[1:20] %>% cbind(coverage, richness, shannon, evenness)
colnames(df_env_div)[21:24] <- c("Coverage", "Richness", "Shannon", "Evenness")
write_tsv(df_env_div, "df_env_div.table", col_names = TRUE)

## Plot rarefaction curve
rarecurve(merge_df[21:ncol(merge_df)], step = 100, xlab = "Number of reads", ylab = "Number of ASVs",
          label = TRUE, xlim=c(0,15000), ylim=c(0,250), lwd=2, lty=1, las=1)

#############################################################
#
#


################# Prepare key dataset for the ASV-based analyses ##################

# "best_hits" table was manually curated to remove ambiguous classification (best hits to ≥2 genus) and to unify the taxonomic rank (genus level)

# Import Edversen's taxonomic assignment information
best_hits <- read_tsv("rep_seqs.blastn.besthit_forR_header.curated.txt", col_names = TRUE)
tax_info <- read_tsv("Edvardsen_hapto_db.tax.txt", col_names = TRUE)

# Merge tax information
asv_tax <- curated_ASVs %>% left_join(best_hits, by = "ASVs") %>% left_join(tax_info, by = "Clone_and_organism")
colnames(asv_tax)[1] <- "name"
asv_tax[is.na(asv_tax)] <- "Unclassified"   # Change NA to "Unclassified"
asv_tax[asv_tax == "unclassified"] <- "Unclassified"    # Change "unclassifies" to "Unclassified" (only for family  category of Pavlomulina)

# Save tax table
write_csv(asv_tax, "asv_tax.csv", col_names = TRUE)

# Make relative abundance table
merge_df2 <- as.data.frame(merge_df)
rownames(merge_df2) <- merge_df2[,1]
hapfreq_df <- merge_df2[,c(-1:-20)]

# Convert to 0-1 ratio
hapfreq_stdz_df <- decostand(hapfreq_df, method="total", MARGIN=1)

#############################################################
#
#


################## PCoA plot with Jacard distance ##################

merge_df <- read_tsv("merge_df.tsv", col_names = TRUE)
ncol_df <- ncol(merge_df)

# Make distance matrix based on the ASV abundance
df_pcoa <- vegdist(merge_df[,21:ncol_df], method="jaccard") 

# Calculate coordinates of NMDS, then extract the coordinates (x and y) of the plot as separated objects.
pc_xy <- cmdscale(df_pcoa, k=2, eig=TRUE)
pc_x <- pc_xy$points[,1]
pc_y <- pc_xy$points[,2]

# Calculate contribution of each axis (eigen value devided by the sum of eigen values)
pc_xy$eig[1:2] / sum(pc_xy$eig)

# Test if the PCoA ordination fits well to the real distance
# https://hoxo-m.hatenablog.com/entry/20120313/p1
eig <- pc_xy$eig           
p1 <- sum(abs(eig[1:2])) / sum(abs(eig)) 
p2 <- sum(eig[1:2] ^ 2) / sum(eig ^ 2)    
cat("Mardia fit measure 1 = ", p1, "\n")
cat("Mardia fit measure 2 = ", p2, "\n")
# Mardia fit measure 1 =  0.5835012
# Mardia fit measure 2 =  0.9464688 <- above the threshold 0.8

# Create df that contains coordinates, env, and diversity estimates of each sample.
metadf <- data.frame(Site=df_env_div[,1],
                     PCoA1 = pc_x,
                     PCoA2 = pc_y,
                     Latitude = df_env_div[,2],
                     Longitude = df_env_div[,3],
                     Temperature = df_env_div[,9],
                     Salinity = df_env_div[10],
                     Chla = df_env_div[13],
                     NO3_NO2 = df_env_div[14],
                     PO4 = df_env_div[15],
                     copy_18S = df_env_div[17],
                     GroupL = df_env_div[19],
                     GroupS = df_env_div[20],
                     Richness = df_env_div[,22],
                     Shannon = df_env_div[,23],
                     Evenness = df_env_div[,24])

# Write table
write_tsv(metadf, "metadf_PCoA.table")

# Draw PCoA plot with Sample group L and richness
nmplot8 <- ggplot(metadf, aes(x=PCoA1, y=PCoA2, fill=GroupL, size=Richness)) +
  geom_point(shape=21) +
  geom_text(aes(x=PCoA1+0.02, label=Site), size=2, hjust=0) +
  scale_color_gradientn(colours = c("blue2", "red2")) +
  theme(legend.position = "right") + theme_bw()
print(nmplot8)
ggsave(path="Figs", filename = "fig.PCoA_GroupL_Richness.pdf", plot= nmplot8, width=125, height=100, units="mm", dpi = 300)

# Add positive and negative of x-axis as a column
metadf_sub <- metadf %>%
  mutate(Pos1 = PCoA1 >= 0) %>% mutate(Pos2 = PCoA2 >= 0)

nmbar1 <- ggplot(metadf_sub, aes(x=Site, y=PCoA1, fill=Pos1)) +
  geom_col(position = "identity", colour= "black", size=0.25) + theme_bw() +
  scale_fill_manual(values = c("#CCEEFF","#FFDDDD"), guide="none") +
  theme(axis.text.x = element_text(angle = -90, hjust = 1))
print(nmbar1)
ggsave(path="Figs", filename = "fig.PCoA_PC1_Site.pdf", plot= nmbar1, width=150, height=50, units="mm", dpi = 300)

#############################################################
#
#


################## Test correlation among variables ##################

# Load data
metadf <- read_tsv("metadf_PCoA.table", col_names = TRUE) #including sample lavel
metadf_select <- as.matrix(metadf[,c(2:11,14:16)])

# Calcurate correlatiuon coefficient and p-value between samples (using the library "Hmisc")
mcor <- rcorr(metadf_select, type = c("spearman"))
# Extract each parameter
r.mcor <- mcor$r
p.mcor <- mcor$P
n.mcor <- mcor$n
# calculate r2 value
r2.mcor <- apply(r.mcor, c(1,2), function(x) {return (x^2)})

# Extract important parameters
r.mcor_select <- as.data.frame(r.mcor[c(5,6,7,8,9),c(1,10,11,12,13)])
r2.mcor_select <- as.data.frame(r2.mcor[c(5,6,7,8,9),c(1,10,11,12,13)])
p.mcor_select <- as.data.frame(p.mcor[c(5,6,7,8,9),c(1,10,11,12,13)])

# Make list for adjusted-q value (Benjamini & Hochberg method)
q.df = NULL
for (i in 1:5){
  q.list <- as.data.frame(p.adjust(as.numeric(p.mcor_select[,i]), method = "BH"))
  colnames(q.list) <- "q"
  q.df <- append(q.df, q.list)
}
# Construct df
q.mcor_select <- data.frame(PCoA1 = q.df[1],
                            copy_18S = q.df[2],
                            Richness = q.df[3],
                            Shannon = q.df[4],
                            Evenness = q.df[5])
# Change row and col names
colnames(q.mcor_select) <- colnames(p.mcor_select)
rownames(q.mcor_select) <- rownames(p.mcor_select)

# Merge r2 and q value
corstat_merged <- data.frame(rownames(r2.mcor_select),
                             r2.mcor_select,
                             q.mcor_select)
colnames(corstat_merged)[1] <- "Expval"

# Save data
write_tsv(r.mcor_select, "r.mcor_select.table")
write_tsv(corstat_merged, "corstat_merged.table")

#############################################################
#
#


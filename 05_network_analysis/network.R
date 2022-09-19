##############KH15-4 haptophyte ASVs analysis for ecological network##############

#Import library
library(magrittr)
library(vegan)
library(gplots)
library(ggplot2)
library(RColorBrewer)
library(MASS)
library(tidyverse)
library(igraph)
library(tidygraph) # For network visualization
library(ggraph) # For network visualization
library(ggforce) # For highlight network cluster
library(graphlayouts) # For highlight network cluster
library(oaqc) # For highlight network cluster
library(pals) # Color palette for many nodes
library(khroma) # For making color pallet
library(ggExtra) # Add density plot to the ggplot
library(Hmisc) # To calculate pairwise r values 'as a matrix"
library(tcltk) # To visualize the progress of "for" function
library(spaa) # For niche analysis

#set directory
print("set working directory")

#load abundance table
merge_df <- as.data.frame(read_tsv("merge_df.tsv", col_names = TRUE))
rownames(merge_df) <- merge_df[,1]
hapfreq_df <- merge_df[,c(-1:-20)]
hap18S_df <- as.data.frame(merge_df[,17])

# Load tax information
asv_tax <- read.table("asv_tax.csv", header=T, sep = ",")
colnames(asv_tax)[1] <- "name"
asv_abb <- read_tsv("asv_abb.txt", col_names = TRUE)

# Load diversity indices
df_env_div <- read_tsv("df_env_div.table", col_names = TRUE)
# Load metadf including PCoA axes
metadf <- read_tsv("metadf_PCoA.table", col_names = TRUE)

# Count number of sample emerged 
sample_emerged_all <- as.numeric(specnumber(hapfreq_df, MARGIN = 2))
num_emerged_all <- data.frame(name=colnames(hapfreq_df), emergence = sample_emerged_all)
hist(sample_emerged_all, breaks=seq(0,46,1))

######################################################
#
#


################# Make absolute abundance profile for m5 and m10 #########################
## As compositional data is not suitable for Pearson and Spearman's metrix, estimating absolute abundance for each ASVs here.
## Reference: Hirano and Yakemoto, BMC Bioinformatics, 2019, doi: 10.1186/s12859-019-2915-1

# Convert to 0-1 ratio
hapfreq_stdz_df <- decostand(hapfreq_df, method="total", MARGIN=1)
write_tsv(hapfreq_stdz_df, "hapfreq_stdz_df.tsv", col_names = TRUE)

# 各サンプル行の要素全てに18S copyを掛け算する（sewwp関数の掛ける値にapplyで18S abundanceを指定。）
hap_abund_df <- sweep(hapfreq_stdz_df, 1, apply(hap18S_df,1, mean), FUN="*")
# Save df
write_tsv(hap_abund_df, "hap_abund_df.tsv", col_names = TRUE)

# Remove ASVs emerged less than 5 samples (for network analysis)
start <- 1
end <- ncol(hap_abund_df)
delete <- NULL
for (i in start:end){
  if (sum(hap_abund_df[,i]!=0) < 5){
    delete <- c(delete, i)
  }
}
hap_abund_m5_df <- hap_abund_df[,-delete]

# Count total reads for each sample
rowSums(hap_abund_m5_df)
# Write data (for network analysis)
write_tsv(hap_abund_m5_df, "hap_abund_m5_df.tsv", col_names = TRUE)
#write_tsv(as.data.frame(colnames(hap_abund_m5_df)), "hap_asv_ids_m5.txt", col_names = FALSE)
# Count number of sample emerged -> used as a size of network nodes
sample_emerged_m5 <- specnumber(hap_abund_m5_df, MARGIN = 2)
num_emerged_m5 <- data.frame(name=colnames(hap_abund_m5_df), emergence = sample_emerged_m5)
barplot(sample_emerged_m5)

# Save gene ids 
m5_ids <- as.data.frame(colnames(hap_abund_m5_df))

##################################################################
#
#


############## Network analysis with Pearson m5 ##################
# Here I use the absolute abundance profile for the Pearson correlation.
# Keep the q values < 0.0001 after BH p-value adjustment.

# Load data
hap_abund_m5_df <- read_tsv("hap_abund_m5_df.tsv")

# Log(x+1) transformation 
hap_abund_m5_df_log <- log1p(hap_abund_m5_df)
rowSums(hap_abund_m5_df_log)

# Convert "charactor" matrix to "numeric" matrix
hap_abund_m5_df_n <- apply(hap_abund_m5_df_log, 2, as.numeric)

# Calculate correlation coefficient and p-value between samples (using the library "Hmisc")
cor.mat.Hmisc <- rcorr(hap_abund_m5_df_n, type = c("pearson"))
r.mat <- cor.mat.Hmisc$r
p.mat <- cor.mat.Hmisc$P
n.mat <- cor.mat.Hmisc$n # should be 0

# Remove upper triangle matrix
r.mat[upper.tri(r.mat, diag = TRUE)] <- NA
p.mat[upper.tri(p.mat, diag = TRUE)] <- NA

# Convert to Long-Format, then merge matrices
# r value
r.list <- r.mat %>%
  as.data.frame() %>%
  mutate(ids1 = colnames(hap_abund_m5_df)) %>%
  gather(key = ids2, value = cor, -ids1)
# p value
p.list <- p.mat %>%
  as.data.frame() %>%
  mutate(ids1 = colnames(hap_abund_m5_df)) %>%
  gather(key = ids2, value = p, -ids1)
# Merge r and p list 
all.list <- cbind(r.list, p.list$p)
colnames(all.list)[4] <- "p"

# FDR (Benjamini-Hochberg) correction for the q value 
q.list <- as.data.frame(p.adjust(as.numeric(p.list[,3]), method = "BH"))
colnames(q.list) <- "q"

# Add q value to the list
all.list.bh <- cbind(all.list, q.list)

# Remove NA and non-significant (q>0.001) nodes
d <- all.list.bh %>% filter(!is.na(cor) & q <= 0.0001)
# Devide into positive and negative
d.posi <- d  %>% filter(cor >= 0)
d.nega <- d  %>% filter(cor <= 0)

# percentage of connection against all possible connections
nrow(d)*100/choose(n=ncol(hap_abund_m5_df), k=2)  # choose: calculate possible sets from n elements by picking k elements up.
nrow(d.posi)*100/choose(n=ncol(hap_abund_m5_df), k=2)
nrow(d.nega)*100/choose(n=ncol(hap_abund_m5_df), k=2)


# Convert df to tbl_graph object
g <- as_tbl_graph(d.posi, directed = FALSE) # if directed graph, chose TRUE

### Calculate network indices with igraph function
# Calculate density
g %>% igraph::graph.density() # -> 除かれた39 ASVがpossible pairsに考慮されていないので値が異なる。
# Add some network features 
g %>% igraph::transitivity()
g %>% igraph::reciprocity()
g <- g %>% mutate(betweenness = centrality_betweenness()) %>% 
  mutate(closeness = centrality_closeness()) %E>% 
  mutate(betweenness = centrality_edge_betweenness())

# Detect community with fast_greedy¥ algorithm
g <- g %N>% mutate(community = as.factor(group_fast_greedy(weights = cor)))

# Evaluate modularity  
gcom <- cluster_fast_greedy(g, weights = E(g)$cor)
dendPlot(gcom, mode="hclust")
plot(gcom, g)
class(gcom)
length(gcom) 
modularity(gcom) 

# Add some taxonomic information 
g <-  g %>% mutate(degree = centrality_degree()) %>% 
  mutate(hub_score = hub_score(g, weights=NA)$vector) %>% 
  left_join(asv_abb, by = "name") %>% 
  left_join(num_emerged_all, by = "name") %>% 
  left_join(asv_tax, by = "name")
# Save as edgelist network file
write.graph(g, "hap_abund_m5_network_pearson.txt", format = "pajek")

# Make table considering Presence/Absence of ASVs for each station and merge it to the g
hapfreq_df_t <- as.data.frame(t(hapfreq_df))
hapfreq_df_t$name <- rownames(hapfreq_df_t)
g_sub <- as.data.frame(g)
g_sub <- g_sub[,c(1,4)]
hapfreq_df_t_select <- g_sub %>% left_join(hapfreq_df_t, by = "name")

# Put community # if the ASV exist, while add "7" if absent
binary_df <- NULL
for (i in 3:ncol(hapfreq_df_t_select)){
  binary = as.data.frame(ifelse(hapfreq_df_t_select[,i] > 0, paste("Module", hapfreq_df_t_select[,2], sep=" "), "Not detected"))
  colnames(binary) <- colnames(hapfreq_df_t_select)[i]
  binary_df <- as.data.frame(append(binary_df, binary))
}
binary_df$name <- hapfreq_df_t_select[,1]
# Add tax information
g <- g %>% left_join(binary_df, by = "name")

# Convert df and save
ggg <- as.data.frame(g)
write_tsv(ggg, "hap_abund_m5_network_pearson2.table", col_names = TRUE)

# Make df specify the color
col_df <- ggg[,17:ncol(ggg)] %>% 
  dplyr::mutate_all(~str_replace(.,pattern="Module 1",replacement = "#C20088")) %>%
  dplyr::mutate_all(~str_replace(.,pattern="Module 2",replacement = "#0075DC")) %>% 
  dplyr::mutate_all(~str_replace(.,pattern="Module 3",replacement = "#993F00")) %>% 
  dplyr::mutate_all(~str_replace(.,pattern="Module 4",replacement = "#4C005C")) %>% 
  dplyr::mutate_all(~str_replace(.,pattern="Module 5",replacement = "#005C31")) %>% 
  dplyr::mutate_all(~str_replace(.,pattern="Module 6",replacement = "#000000")) %>% 
  dplyr::mutate_all(~str_replace(.,pattern="Not detected",replacement = "#FFFFFD"))


## Define node color: the custom function using Color Brewer
colnodes <- as.vector(alphabet(26)) # in package "pals"
colnodes.pearson <- c("#C20088","#0075DC", "#993F00","#4C005C","#005C31", "#000000")
colnode_circle <- c("gray14","gray14","gray14","gray14","gray14","gray14","gray14","gray14","gray14","gray14","gray14",
                    "gray14","gray14","gray14","gray14","gray14","gray14","gray14","gray14","gray14","gray14","gray14",
                    "gray14","gray14","gray14","gray14","gray14","gray14","gray14","gray14","gray14","gray14","gray14")

### Visualize network
# Random
plot1 <- g %>%
  ggraph(layout = "kk") +  # chose kk, fr, or graphopt
  geom_edge_link(aes(width=cor), alpha = 0.6, colour = "gray50") +
  scale_edge_width(range = c(0.1, 1)) +
  scale_color_manual(values = colnode_circle) +
  scale_fill_manual(values = colnodes.pearson) +
  geom_node_point(aes(size = degree, colour = community, fill=community), shape = 21, alpha=0.7) +
  #geom_node_text(aes(label = name, colour=community,), repel = TRUE) +
  theme_graph(background = "white")
plot1
#保存
quartz(type="pdf", width=13, height=10, dpi = 300, file = "Fig.hap_abund_m5_pearson.e-4_network_kk_link.pdf")
plot1
dev.off()

##################################################################
#
#


########### Correlation between phylogenetic distance and co-abundance trend m5 ###############

## Here I compare co-abundance trend (pearson correlation) and phylogenetic distance for selected OTUs.

## 1 ## Prepare genetic distance
# Phylogenetic distance (%) can be calculated by "Geneious" function from alignment file
## Then convert matrix to list
dist.mat <- read_tsv("m5_dist.txt", col_names = TRUE)
rownames(dist.mat) <- dist.mat$X1
dist.mat <- as.matrix(dist.mat[,-1])
dist.mat[upper.tri(dist.mat, diag = TRUE)] <- NA
gendist_all <- dist.mat %>%
  as.data.frame() %>%
  mutate(ids1 = colnames(hap_abund_m5_df)) %>%
  gather(key = ids2, value = gendist, -ids1)

## 2 ## Merge with the Pearson's correlation data
# Here we use "all.list" generated in the previous section.
## Merge both matrix
dist_cor_merged <- cbind(gendist_all, all.list)
dist_cor_merged2 <- dist_cor_merged[,c(-4,-5,-7)]

# Save df
write_tsv(dist_cor_merged2, "dist_cor_merged2.table", col_names = TRUE)


# Add taxonomic info
## Make sub tax 
asv_tax_sub1 <- asv_tax[, c(1,6:8)]
colnames(asv_tax_sub1) <- c("ids1", "Family1", "Genus1", "Species1")
asv_tax_sub2 <- asv_tax[, c(1,6:8)]
colnames(asv_tax_sub2) <- c("ids2", "Family2", "Genus2", "Species2")

# Merge tax info for both nodes and filter "NA" correlation out
dist_cor_merged3 <- dist_cor_merged2 %>% 
  left_join(asv_tax_sub1, by="ids1") %>% 
  left_join(asv_tax_sub2, by="ids2") %>% 
  filter(!is.na(gendist))        # Remove rows with gendist == NA


# Specify genus which will be parsed later
# selected 11 genera
genus_list <- c("Chrysochromulina","Phaeocystis","Clade_D","Prymnesium","Clade_B4","Clade_HAP4",
                "Emiliania_Gephyrocapsa","Syracosphaera","Clade_HAP5","Clade_HAP2","Pavlomulina")

# List of 30 genera 
genus_list_all_remuc <- c("Chrysochromulina","Phaeocystis","Clade_D","Prymnesium","Clade_B4","Clade_HAP4","Emiliania_Gephyrocapsa",
                    "Syracosphaera","Clade_HAP5","Clade_HAP2","Pavlomulina","Clade_E","Dicrateria","Algirosphaera","Haptolina",
                    "Braarudosphaera","Clade_HAP3","Clade_B3","Chrysoculter","Scyphosphaera","Clade_F","Umbilicosphaera","Chrysotila",
                    "Gladiolithus","Calcidiscus","Oolithus","Exanthemachrysis","Clade_B5","Coronosphaera","Pavlova")


# !!!!!!!! This process requires huge time !!!!!!!!!

# Set the progress visualization tool for "for" function
pb <- tcltk::txtProgressBar(min = 1, max = nrow(dist_cor_merged3), style = 3)
## Classify edges with and without connecting to the specific genus
# Make cols to classify if each node is derived from major genara 
Genus_all <- NULL 
for (j in genus_list_all_remuc){
  genus <- data.frame() 
  for (i in 1:nrow(dist_cor_merged3)){
    setTxtProgressBar(pb, i) # Visualize progress for each of the search
    if (dist_cor_merged3[i,6] == j && dist_cor_merged3[i,9] == j){ 
      genus[i,1] <- "Intra-genus"
    } else if (dist_cor_merged3[i,6] == j || dist_cor_merged3[i,9] == j){ 
      genus[i,1] <- "Inter-genus"
    } else{
      genus[i,1] <- "Other"    
    }
  }
  Genus_all <- append(Genus_all, genus)
}
dist_cor_merged4 <- cbind(dist_cor_merged3, Genus_all)
colnames(dist_cor_merged4)[11:40] <- genus_list_all_remuc
colnames(dist_cor_merged4)[4] <- "Pearson"

# Make category column for all taxa
Edge_all <- data.frame()
for (i in 1:nrow(dist_cor_merged4)){
  if ("Intra-genus" %in% dist_cor_merged4[i, c(11:40)]){
    Edge_all[i,1] <- "Intra-genus"
  } else if ("Inter-genus" %in% dist_cor_merged4[i, c(11:40)]){
    Edge_all[i,1] <- "Inter-genus"
  } else {
    Edge_all[i,1] <- "Other"
  }
}
dist_cor_merged4 <- cbind(dist_cor_merged4, Edge_all)
colnames(dist_cor_merged4)[41] <- "Total"

# Write df
#write_tsv(dist_cor_merged4, "dist_cor_merged4_m5.table", col_names = TRUE)

# Load df
#dist_cor_merged4 <- read_tsv("dist_cor_merged4_m5.table", col_names = TRUE)


#### Plot relationship between genetic distance and network correlation

# Define color to classify tax
col_tax  <- c("indianred", "forestgreen", "cornflowerblue")

# For all with color
# Remove ambiguous connection
dist_total_select <- dist_cor_merged4 %>% filter(dist_cor_merged4$Total != "Other") # Filter out unrelated interaction
m1 <- ggplot(dist_total_select, aes(x=gendist, y=Pearson, color=Total)) + theme_bw() +
  ylim(-1,1) +
  geom_point(size=1, shape=16, alpha=0.03) + 
  scale_color_manual(breaks = c("Intra-genus", "Inter-genus"), values = col_tax) +
  theme(legend.position="bottom", legend.direction = "horizontal") +
  stat_smooth(method=loess, size = 0.5, alpha=0.3, fullrange = FALSE, se = FALSE)
# Add density plot
mm1 <- ggMarginal(m1, type = "density", size = 8, margins = "both", groupColour = TRUE, groupFill = TRUE) 
print(mm1)
ggsave(path="Figs_m5", filename = "Fig.gendist_Pearson_total_density_colored.pdf", plot= mm1, width=110, height=125, units="mm", dpi = 300)


# Statistical test for each plot
df_same_genus <- dist_cor_merged4 %>% filter(Total == "Intra-genus") #7528 obs
df_diff_genus <- dist_cor_merged4 %>% filter(Total == "Inter-genus") #42842 obs

cor.test(df_same_genus$gendist, df_same_genus$Pearson, method = "pearson")
# cor 0.1102565, t = 9.6237, df = 7526, p-value < 2.2e-16 
cor.test(df_diff_genus$gendist, df_diff_genus$Pearson, method = "pearson")
# cor 0.01520026, t = 3.1465, df = 42840, p-value = 0.001654

### Select test to be used 
# Shapiro-Wilk test to check if two variable is normal distributed
shapiro.test(x=df_same_genus$Pearson[1:5000])
shapiro.test(x=df_diff_genus$Pearson[1:5000])
#W = 0.98881, p-value < 2.2e-16
#W = 0.9878, p-value < 2.2e-16
# -> Non-parametric method should be used (Mann-Whitney U test)

# F-test to check if two variable is homoscedasticity distributed (等分散)
var.test(x=df_same_genus$Pearson, y=df_diff_genus$Pearson)
# F = 1.0237, num df = 7527, denom df = 42841, p-value = 0.1833
t.test(x=df_same_genus$Pearson, y=df_diff_genus$Pearson, var.equal=F,paired=F)

# Mann-Whitney U test
wilcox.test(df_same_genus$Pearson, df_diff_genus$Pearson, correct=FALSE)
# W = 185677480, p-value < 2.2e-16

# Welch's t-test
t.test(df_same_genus$Pearson, df_diff_genus$Pearson, var.equal=F,paired=F)
# t = 21.428, df = 10280, p-value < 2.2e-16

##################################################################
#
#


########## Genus-level statistics m5 #######

# Generate a table summaryzing the following info
#1. Count number of ASVs in the curated table for each genus. 
#2. Count nunber of ASVs in the network for each genus.
#3. Calcurate mean contribution (%) of each genus across samples.

#Load df of curated data
asv_tax <- read.table("asv_tax.csv", header=T, sep = ",")
colnames(asv_tax)[1] <- "name"

# Load df of network features
gg_net <- read_tsv("hap_abund_m5_network_pearson.table")

# Load genus contribution 
genus_merged_stdz <- read_tsv("asv_curate_stdz_genus_merged.table")

# List of 30 genera + 1 Unclassified
genus_list_all <- c("Chrysochromulina","Phaeocystis","Clade_D","Prymnesium","Clade_B4","Clade_HAP4","Emiliania_Gephyrocapsa",
                    "Syracosphaera","Clade_HAP5","Clade_HAP2","Pavlomulina","Clade_E","Dicrateria","Algirosphaera","Haptolina",
                    "Braarudosphaera","Clade_HAP3","Clade_B3","Chrysoculter","Scyphosphaera","Clade_F","Umbilicosphaera","Chrysotila",
                    "Gladiolithus","Calcidiscus","Oolithus","Exanthemachrysis","Clade_B5","Coronosphaera","Pavlova","Unclassified")

# Execute command
genus_summary <- NULL
for (i in genus_list_all){
  ASVall <- sum(str_count(asv_tax$Genus, i))
  ASVnet <- sum(str_count(gg_net$Genus, i))
  Net_ratio <- ASVnet*100/ASVall
  rownum <- which(genus_merged_stdz[,] == i) # Specify row number
  Meancont <- mean(as.numeric(genus_merged_stdz[rownum, 2:46]))
  Maxcont <- max(as.numeric(genus_merged_stdz[rownum, 2:46]))
  Mincont <- min(as.numeric(genus_merged_stdz[rownum, 2:46]))
  df <- cbind(i, ASVall, ASVnet, Net_ratio, Meancont, Maxcont, Mincont)
  genus_summary <- rbind(genus_summary, df)
}
genus_summary <- as.data.frame(genus_summary, stringsAsFactors = FALSE)
colnames(genus_summary) <- c("Genus", "ASVall", "ASVnet", "Net_ratio", "Meancont", "Maxcont", "Mincont")

# Adjust type of df as these are stored as characters
genus_summary <- data.frame(Genus = as.character(genus_summary$Genus),
                            ASVall = as.numeric(genus_summary$ASVall),
                            ASVnet = as.numeric(genus_summary$ASVnet),
                            Net_ratio = as.numeric(genus_summary$Net_ratio),
                            Meancont = as.numeric(genus_summary$Meancont),
                            Maxcont = as.numeric(genus_summary$Maxcont),
                            Mincont = as.numeric(genus_summary$Mincont)
)

# Save data
write_tsv(genus_summary, "genus_summary.table", col_names = TRUE)

# Load data
genus_summary <- read_tsv("genus_summary.table")

# Plot 
ballp1 <- ggplot(genus_summary, aes(x=Meancont, y=ASVall, size=ASVnet)) +
  geom_point(shape=21, colour="black", fill="cornsilk") + 
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  scale_size_area(max_size=8) + 
  annotation_logticks() +  # Add logarithmic scale
  scale_size(range = c(1, 8)) +
  scale_x_log10() +
  scale_y_log10() +
  ggrepel::geom_text_repel(ggplot2::aes(label = Genus), size = 2.5)
print(ballp1)
ggsave(path="Figs_m5", filename = "Fig.meancont_richness_paper.pdf", plot= ballp1, width=150, height=125, units="mm", dpi = 300)
ggsave(path="Figs_m5", filename = "Fig.meancont_richness_paper_L.pdf", plot= ballp1, width=123, height=100, units="mm", dpi = 300)

##################################################################
#
#


########## ASV-level niche and network feature analysis m5 #######

# Load ASV contribution with the values 0-1 (without rownames)
hapfreq_stdz_df <- read_tsv("hapfreq_stdz_df.tsv")
hap_abund_m5_df <- read_tsv("hap_abund_df.tsv")   # Used for Levins niche width calculation

#Load df of curated data
asv_tax <- read.table("asv_tax.csv", header=T, sep = ",")
colnames(asv_tax)[1] <- "name"

# Load df of network features
gg_net <- read_tsv("hap_abund_m5_network_pearson.table")

# Calcurate Levins' niche breadth
ASV_niche <- niche.width(hap_abund_m5_df, method = "levins")
ASV_niche <- as.data.frame(t(ASV_niche))
# Standardize Levins' index to 0-1 value 
ASV_niche_stdz <- as.data.frame(apply(ASV_niche, c(1,2), function(x) {print((x-1)/(46-1))}))
#Merge
ASV_niche <- data.frame(as.data.frame(rownames(ASV_niche)), ASV_niche_stdz) 
colnames(ASV_niche) <- c("name", "Levins")
# Count mean relative abundance (%) for each ASVs, then add vector as a column
ASV_niche[, "Meancont"] <- as.data.frame(apply(hapfreq_stdz_df[, 1:437]*100, 2, mean))
ASV_niche[, "Maxcont"] <- as.data.frame(apply(hapfreq_stdz_df[, 1:437]*100, 2, max))
ASV_niche[, "Mincont"] <- as.data.frame(apply(hapfreq_stdz_df[, 1:437]*100, 2, min))

# Make merged df for ASVs
ASV_summary <- ASV_niche %>% left_join(num_emerged_all, by="name") %>% 
               left_join(asv_tax, by="name") %>% 
               left_join(gg_net[,1:7], by="name")

# List of selected genus
genus_list <- c("Chrysochromulina","Phaeocystis","Clade_D","Prymnesium","Clade_B4","Clade_HAP4",
                "Emiliania_Gephyrocapsa","Syracosphaera","Clade_HAP5","Clade_HAP2","Pavlomulina")
# Make column with selected genus names
Genus_select <- data.frame(Genus_select = rep("Other", 437)) # Make empty df
Genus_select$Genus_select <- as.character(Genus_select$Genus_select) # R6.**だとdfの型がfactorになってしまうので変換
for (i in genus_list){
  for (j in 1:437){
    if (ASV_summary[j,12] == i){
      Genus_select[j,] <- i
    }
  }
}
# Add as a new column
ASV_summary <- data.frame(ASV_summary, Genus_select=Genus_select)

# Make column defining ecological strategy based on Levins niche width
# The definition is original, although Logares et al(2013) used much less values for generalist
# 上位1/3などで定義する?
strategy <- data.frame(Strategy = rep(NA, 437)) # Make empty df
for (i in 1:437){
  if (ASV_summary[i,2] > 0.6){
    strategy[i,] <- "Generalist"
  } else if (ASV_summary[i,2] < 0.4){
    strategy[i,] <- "Specialist"
  } else {
    strategy[i,] <- "intermediate"
  }
}
# Add as a new column
ASV_summary <- data.frame(ASV_summary, Strategy=strategy)


# Save data
write_tsv(ASV_summary, "ASV_summary.table", col_names = TRUE)
# Read table (Can be restart from here)
ASV_summary <- read_tsv("ASV_summary.table")


# Max contrbution vs Levins nitch width for paper
cor.test(ASV_summary$Maxcont, ASV_summary$Levins, method = "spearman")
ballp2 <- ggplot(ASV_summary, aes(x=Maxcont, y=Levins)) +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  geom_point(shape=20, colour="black", alpha=0.8) + 
  scale_size_area(max_size=3) +
  annotation_logticks(sides = "b") +  # Add logarithmic scale, only x axis
  scale_x_log10(limits= c(0.01,100)) 
print(ballp2)
ggsave(path="Figs_m5", filename = "Fig.maxcont_Levins_paper_L.pdf", plot= ballp2, width=98, height=95, units="mm", dpi = 300)

## Levins' standardized index across genera
boxp2 <- ggplot(ASV_summary, aes(x=Genus_select, y=Levins, fill=Genus_select)) +
  geom_boxplot() + theme_bw() +
  scale_y_continuous(limits = c(0, 0.8)) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0)) +
  geom_jitter(shape=16, position=position_jitter(0.3), size =0.1, alpha=0.1)
print(boxp2)
ggsave(path="Figs_m5", filename = "Fig.gendist_boxplot_all.pdf", plot= boxp2, width=90, height=110, units="mm", dpi = 300)

##################################################################
#
#



############ Network module-based analysis m5 #####################

# Make a plot showing ralative contributions of each community

# Load dfs
ASV_summary <- read_tsv("ASV_summary.table", col_names = TRUE)
ASV_summary_select <- ASV_summary[,c(1,17)]
hapfreq_stdz100_dft　<- read.csv("hapfreq_stdz100_dft_rown.csv", header=T)
colnames(hapfreq_stdz100_dft)[1] <- "name"

# Add community information
hapfreq_stdz100_com <- hapfreq_stdz100_dft %>%  left_join(ASV_summary_select, by = "name")

# Extract freq and genus
hapfreq_stdz100_com2 <- hapfreq_stdz100_com[,c(2:48)]
# grouping by genus classification and summarize them
hapfreq_stdz100_com_merged <- hapfreq_stdz100_com2 %>% group_by(community) %>% summarise_each(sum)
# Further merge communities 4-NA to "Others"
hapfreq_stdz100_com_merged2 <- cbind(hapfreq_stdz100_com_merged, c("Module 1", "Module 2", "Module 3",
                                                                   "Others", "Others", "Others", "Others"))
colnames(hapfreq_stdz100_com_merged2)[48] <- "Module2"
hapfreq_stdz100_com_merged3 <- hapfreq_stdz100_com_merged2[,2:48] %>% group_by(Module2) %>% summarise_each(sum)
hapfrew_com_merged <- as.data.frame(t(hapfreq_stdz100_com_merged3[,2:47]))
hapfrew_com_merged <- cbind(merge_df$Site, hapfrew_com_merged)
colnames(hapfrew_com_merged) <- c("Station","Module 1", "Module 2", "Module 3", "Others")

# Convert to the tydi format
comm_tydi <- gather(hapfrew_com_merged, key = Module, value = Contribution, -Station)
# Define order
comm_tydi$Module <- factor(comm_tydi$Module, levels=c("Others", "Module 3", "Module 2", "Module 1"))
# Define color
col_module2 <- c("gray60", "#b7784c", "#108fff", "#DD6677")

# Plot
compbar_rel <- ggplot(comm_tydi, aes(x=Station, y=Contribution, fill=Module)) +
  geom_col(colour="black", width = 1) + 
  scale_fill_manual(values = col_module2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -90, hjust = 1))
print(compbar_rel)
ggsave(path="Figs_m5", filename = "Fig.com_all_rel.pdf", plot= compbar_rel, width=190, height=50, units="mm", dpi = 300)

### Parse "gg" to obtain count table of each community for each samples
# Load df
gg <-read_tsv("hap_abund_m5_network_pearson.table", col_names = TRUE)

# Make station list for the loop
station <- colnames(gg)[17:62]
# Make count df
com_st_count <- NULL
for (i in 1:6){
  count_sum <- NULL
  for (j in 17:62){
    count <- sum(with(gg, gg$community == i & gg[j] == "Present"))
    count_sum <- append(count_sum, count)
  }
  com_st_count <- as.data.frame(cbind(com_st_count, count_sum))
}
com_st_count <- data.frame(Station = station, com_st_count)
colnames(com_st_count) <- c("Station","Module 1","Module 2","Module 3","Module 4","Module 5","Module 6")
com_st_count_plus <- cbind(com_st_count, merge_df[,c(2,3,17,18,19,20)])

# Make table summarizing the community 
Others_div <- df_env_div[,22] - rowSums(com_st_count[,2:4]) #ASVs not included in the network
com_st_count2 <- data.frame(com_st_count[,1:4], Others_div)
colnames(com_st_count2) <- c("Station", "Module 1", "Module 2", "Module 3", "Others")

# Convert to the tydi format
comm_count_tydi <- gather(com_st_count2, key = Module, value = Richness, -Station)
# Define order
comm_count_tydi$Module <- factor(comm_count_tydi$Module, levels=c("Others", "Module 3", "Module 2", "Module 1"))

# Plot
compbar_abs <- ggplot(comm_count_tydi, aes(x=Station, y=Richness, fill=Module)) +
  geom_col(colour="black", width = 1) + 
  scale_fill_manual(values = col_module2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -90, hjust = 1))
print(compbar_abs)
ggsave(path="Figs_m5", filename = "Fig.com_all_abs.pdf", plot= compbar_abs, width=190, height=50, units="mm", dpi = 300)


### Make plot for PCoA1 vs Modules2+3

# Generate df having PCoA value
st_merge_df <- cbind(hapfrew_com_merged, metadf[,2:5])
st_merge_df <- st_merge_df %>% mutate(ECSmodules = apply(st_merge_df[,3:4], 1, sum))
# Change col names
colnames(st_merge_df)[2:4] <- c("Module_1", "Module_2", "Module_3")

#Plot
pm2 <- ggplot(st_merge_df, aes(x=PCoA1, y=ECSmodules)) +
  ylim(0,100) +
  geom_point(size=2, shape=16) + 
  stat_smooth(method=lm, size = 0.5, se = FALSE, fullrange = FALSE) +
  theme_bw()
print(pm2)
ggsave(path="Figs_m5", filename = "Fig.PCoA1_ECSmodule_bw.pdf", plot= pm2, width=70, height=70, units="mm", dpi = 300)

##################################################################
#
#


###### This script aims to extract environmental data from metadata table ######

# Load library
library(tidyverse)
library(geosphere)

# Set directory
print("set working directory")

# Site data
sites <- read_tsv("sample_sites.tsv", col_names = TRUE) # Site info of DNA samples
amembo_all <- read_tsv("AMEMBO_merged_select.tsv", col_names = TRUE)
ume_env <- read_tsv("umezawa_env.tsv", col_names = TRUE)


## Search closest AMEMBO data by long and lat, and then combine it with st data
st_phys <- data.frame()  
for (i in 1:nrow(sites)){
  site_loc <- as.numeric(c(sites[i,3], sites[i,2]))
  dist_list <- NULL
  for (j in 1:nrow(amembo_all)){
    env_loc <- as.numeric(c(amembo_all[j,2], amembo_all[j,1]))
    dist <- distm(site_loc, env_loc, fun=distGeo)
    dist_list <- append(dist_list, dist)
  }
  min_dist <- min(dist_list)
  min_row <- which.min(dist_list)
  each_info <- cbind(sites[i,1:3], min_row, min_dist, amembo_all[min_row,1:6])
  st_phys <- rbind(st_phys, each_info)
}
st_phys_df <- as.data.frame(st_phys)


## Search closest nutrient data by long and lat, and then combine it with st data
st_phys_chem <- data.frame()  
for (i in 1:nrow(st_phys_df)){
  site_loc <- as.numeric(c(st_phys_df[i,3], st_phys_df[i,2]))
  dist_list <- NULL
  for (j in 1:2465){
    env_loc <- as.numeric(c(ume_env[j,2], ume_env[j,1]))
    dist <- distm(site_loc, env_loc, fun=distGeo)
    dist_list <- append(dist_list, dist)
  }
  min_dist <- min(dist_list)
  min_row <- which.min(dist_list)
  each_info <- cbind(st_phys_df[i,1:11], min_row, min_dist, ume_env[min_row,1:11])
  st_phys_chem <- rbind(st_phys_chem, each_info)
}
dist_result <- as.data.frame(st_phys_chem)
write_tsv(st_phys_chem, "Results/station_phys_chem.table", col_names = TRUE)



# Note that the original gps coordinates are not collect as some of the stations are wrong recorded as follows:
# 56.060 -> 56.60

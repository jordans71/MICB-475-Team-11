#### Install Packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")

library(phyloseq)
library(tidyverse)
library(ggplot2)
library(ape)
library(vegan)

# Load Rdata
load("parkinsons_final_sleep.RData")

# Extracting OTU data
otu_table <- data.frame(t(otu_table(parkinsons_final_sleep)))
otu_table$ID <- rownames(otu_table)

# Extracting metadata
metadata <- data.frame(sample_data(parkinsons_final_sleep))
metadata$ID <- rownames(metadata)

# Load the raw taxonomy file
tax <- data.frame(tax_table(parkinsons_final_sleep))

# Formatting the taxa dataframe and cleaning names
tax_mat <- tax[,-1]
tax_mat <- data.frame(tax_mat)
tax_mat$Phylum <- gsub("^...","",tax_mat$Phylum)
tax_mat$Class <- gsub("^...","",tax_mat$Class)
tax_mat$Order <- gsub("^...","",tax_mat$Order)
tax_mat$Family <- gsub("^...","",tax_mat$Family)
tax_mat$Genus <- gsub("^...","",tax_mat$Genus)
tax_mat$Species <- gsub("^...","",tax_mat$Species)
tax_mat$ASV <- rownames(tax_mat)

# Joining OTU and metadata and taxanomic information
otu_meta <- inner_join(metadata, otu_table, by = "ID")

#transforming the OTU matrix to a single column called abundance. 53 represents the number of metadata columns we want to exclude
?gather
grouped = gather(otu_meta, key = "ASV", value = "abundance", -(1:27))

#Joining the taxa information with otu_meta
grouped_taxa = inner_join(tax_mat, grouped, by = "ASV", multiple = "all")


grouped_taxa$legend = paste(grouped_taxa$sheetwashfreq_binned) #Can add gender for when gender becomes a consideration

#Calculating the relative abundance for each indivudal
ppl <- unique(grouped_taxa$ID) 
data_rel = data.frame() #Create empty dataframe
for (i in ppl){
  df = grouped_taxa %>%
    filter(ID == i)
  
  df_sum = df %>%
    group_by(ID,legend,sex, sheetwashfreq_binned, Phylum) %>%
    summarize(abs = sum(abundance))
  
  df_sum$abs = as.numeric(as.character(df_sum$abs))
  count = sum(df_sum$abs)
  df_sum$rel_abs = df_sum$abs/count*100
  
  data_rel = rbind(data_rel, df_sum)
  
}

data_rel$sheetwashfreq_binned = factor(data_rel$sheetwashfreq_binned, levels = c("low","medium","high")) #create the order for low, medium, high in the plot
ggplot(data =data_rel, aes(ID,rel_abs, fill = Phylum))+ #Generating the plot with X axis equal to individual
  geom_col()+
  theme(axis.text.x = element_text(angle = -90))+
  facet_grid(cols = vars(sheetwashfreq_binned), scales = "free_x", space = "free_x")



#Generating the relative abundance for all individuals within low/medium/high sheetwashing freq.
levels <- unique(grouped_taxa$sheetwashfreq_binned)
data_rel = data.frame()
for (i in levels){
  
  df = grouped_taxa %>%
    filter(sheetwashfreq_binned == i)
  
  df_sum = df %>%
    group_by(ID,legend,sex, sheetwashfreq_binned, Phylum) %>%
    summarize(abs = sum(abundance))
  
  df_sum$abs = as.numeric(as.character(df_sum$abs))
  count = sum(df_sum$abs)
  df_sum$rel_abs = df_sum$abs/count*100
  
  data_rel = rbind(data_rel, df_sum)
  
}

#This plot represents the average relative abundance for each phylum across the different sheetwash frequency levels.
data_rel$sheetwashfreq_binned = factor(data_rel$sheetwashfreq_binned, levels = c("low","medium","high")) #create the order for low, medium, high in the plot
ggplot(data =data_rel, aes(sheetwashfreq_binned,rel_abs, fill = Phylum))+#Generating the plot with X axis equal to sheetwash_freq_binned
  geom_col()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = -90),
        legend.position = "none")


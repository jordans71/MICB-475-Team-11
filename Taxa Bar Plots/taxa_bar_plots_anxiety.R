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
load("parkinsons_final_anxiety.RData")

# Extracting OTU data
otu_table <- data.frame(t(otu_table(parkinsons_final_anxiety)))
otu_table$ID <- rownames(otu_table)

# Extracting metadata
metadata <- data.frame(sample_data(parkinsons_final_anxiety))
metadata$ID <- rownames(metadata)

# Load the raw taxonomy file
tax <- data.frame(tax_table(parkinsons_final_anxiety))

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

#transforming the OTU matrix to a single column called abundance. 
grouped = gather(otu_meta, key = "ASV", value = "abundance", -(1:102))

#Joining the taxa information with otu_meta
grouped_taxa = inner_join(tax_mat, grouped, by = "ASV", multiple = "all")


grouped_taxa$legend = paste(grouped_taxa$Disease,grouped_taxa$anxitey_binned) 


#Generating the relative abundance.
levels <- unique(grouped_taxa$legend)
data_rel = data.frame()
for (i in levels){
  
  df = grouped_taxa %>%
    filter(legend == i)
  
  df_sum = df %>%
    group_by(anxitey_binned,legend, Disease, Phylum) %>%
    summarize(abs = sum(abundance))
  
  df_sum$abs = as.numeric(as.character(df_sum$abs))
  count = sum(df_sum$abs)
  df_sum$rel_abs = df_sum$abs/count*100
  
  data_rel = rbind(data_rel, df_sum)
  
}

data_rel_sum = data_rel %>%
  group_by(legend,Disease,anxitey_binned,Phylum) %>%
  summarise(mean_rel_abs = sum(rel_abs))

#Filtering for phyla that represent relative abundance greater than 1% of the group.
data_rel_sum_filtered = data_rel_sum %>%
  filter(mean_rel_abs> 1)

#This plot represents the average relative abundance for each phylum across different levels.
data_rel_sum_filtered$Disease = factor(data_rel_sum_filtered$Disease, levels = c("Control","PD")) #create the order for control, PD in the plot
ggplot(data =data_rel_sum_filtered, aes(anxitey_binned,mean_rel_abs, fill = Phylum))+
  geom_col()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = -90),
        axis.text = element_text(size = 15, face = "bold"),
        strip.text = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(size = 15, face = "bold")
        ,legend.text = element_text(size = 10, face = "bold"),
        legend.title = element_text(size = 15, face = "bold"))+
  facet_grid(cols = vars(Disease), scales = "free_x", space = "free_x")+
  labs(x = "", y = "Relative abundance (%)")

#####################################################################
#CODE FOR PLOTS ABOUT SPECIFIC PHYLUM
####################################################################

#Generating relative abundance at the genus level

#Define the control/PD + yes/no levels
levels <- unique(grouped_taxa$legend)
#Create a new empty dataframe.
data_rel_genus = data.frame()

#Run a loop to generate the relative abundance of each genus for different groups.
for (i in levels){
  
  df = grouped_taxa %>%
    filter(legend == i)
  
  df_sum = df %>%
    group_by(anxitey_binned,legend,Disease,Phylum,Order,Family,Class, Genus) %>%
    summarize(abs = sum(abundance))
  
  df_sum$abs = as.numeric(as.character(df_sum$abs))
  count = sum(df_sum$abs)
  df_sum$rel_abs = df_sum$abs/count*100
  
  data_rel_genus = rbind(data_rel_genus, df_sum)
  
}

#Filter for only Bacteroidota 
data_rel_proteobacteria = data_rel_genus %>%
  filter(Phylum == "Bacteroidota") %>% #Keep only bacteria that are in the Proteobacteria phylum. Change this line for different phyla.
  group_by(legend,Disease,anxitey_binned,Phylum, Order, Class, Family,Genus) %>%
  summarise(mean_rel_abs = sum(rel_abs))%>% #Add relative abundances together for multiple species that have the same genus
  filter(mean_rel_abs>= 1) #Remove Genus that have a relative abundance less than 1%

data_rel_proteobacteria$Disease = factor(data_rel_proteobacteria$Disease, levels = c("Control","PD")) #create the different levels in the plot
ggplot(data =data_rel_proteobacteria, aes(anxitey_binned,mean_rel_abs, fill = Genus))+
  geom_col()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = -90),
        axis.text = element_text(size = 15, face = "bold"),
        strip.text = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(size = 15, face = "bold")
        ,legend.text = element_text(size = 10, face = "bold"),
        legend.title = element_text(size = 15, face = "bold"))+
  facet_grid(cols = vars(Disease), scales = "free_x", space = "free_x")+
  labs(x = "", y = "Relative abundance (%)")

#Filter for only Firmicutes 
data_rel_proteobacteria = data_rel_genus %>%
  filter(Phylum == "Firmicutes") %>% #Keep only bacteria that are in the Proteobacteria phylum. Change this line for different phyla.
  group_by(legend,Disease,anxitey_binned,Phylum, Order, Class, Family,Genus) %>%
  summarise(mean_rel_abs = sum(rel_abs))%>% #Add relative abundances together for multiple species that have the same genus
  filter(mean_rel_abs>= 1) #Remove Genus that have a relative abundance less than 1%

data_rel_proteobacteria$Disease = factor(data_rel_proteobacteria$Disease, levels = c("Control","PD")) #create the different levels in the plot
ggplot(data =data_rel_proteobacteria, aes(anxitey_binned,mean_rel_abs, fill = Genus))+
  geom_col()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = -90),
        axis.text = element_text(size = 15, face = "bold"),
        strip.text = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(size = 15, face = "bold")
        ,legend.text = element_text(size = 10, face = "bold"),
        legend.title = element_text(size = 15, face = "bold"))+
  facet_grid(cols = vars(Disease), scales = "free_x", space = "free_x")+
  labs(x = "", y = "Relative abundance (%)")


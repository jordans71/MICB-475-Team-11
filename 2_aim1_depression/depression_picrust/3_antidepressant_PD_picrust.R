library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)
library(DESeq2)

source("DESeq2_function.R")


#Importing the pahtway PICrsut2
abundance_file <- "Picrust analysis _path_abun_unstrat.tsv"
abundance_data <- read_delim(abundance_file, delim = "\t", col_names = TRUE, trim_ws = TRUE)
abundance_data  =as.data.frame(abundance_data)


#rownames(abundance_data_1) = abundance_data_1$pathway
#abundance_data = abundance_data_1[,-1]

metadata <- read_delim("Picrust analysis _parkinsons_metadata_new_edited.csv")

library("ggh4x")

#Example Looking at depression within PD patients

#Filter your metadata as needed to look at specific comparisons
PD_metadata = metadata %>%
  filter(Disease == "PD")

#Remove NAs for anxitey
PD_metadata = PD_metadata[!is.na(PD_metadata$Antidepressant_use),]

#Filtering the abundance table to only include samples that are in the filtered metadata
sample_names = PD_metadata$`X.SampleID`
sample_names = append(sample_names, "pathway")
abundance_data_filtered = abundance_data[, colnames(abundance_data) %in% sample_names] #This step is the actual filtering

#Removing individuals with no data that caused a problem for pathways_daa()
abundance_data_filtered =  abundance_data_filtered[, colSums(abundance_data_filtered != 0) > 0]

#Filtering out columns that represent a total abundance < 10000
#Filtering out rows (pathways) that have a total count < 100
abundance_data_filtered = abundance_data_filtered[,colSums(abundance_data_filtered[,-1]) > 10000]
abundance_data_filtered = abundance_data_filtered[rowSums(abundance_data_filtered[,-1]) > 100,]

#Ensuring the rownames for the abundance_data_filtered is empty. This is required for their functions to run.
rownames(abundance_data_filtered) = NULL

#verify samples in metadata match samples in abundance_data
abun_samples = rownames(t(abundance_data_filtered[,-1])) #Getting a list of the sample names in the newly filtered abundance data
PD_metadata = PD_metadata[PD_metadata$`X.SampleID` %in% abun_samples,] #making sure the filtered metadata only includes these samples


#Perform pathway DAA using LinDA method
abundance_daa_results_df <- pathway_daa(abundance = abundance_data_filtered %>% column_to_rownames("pathway"), metadata = PD_metadata, group = "Antidepressant_use", daa_method = "DESeq2")

# Annotate MetaCyc pathway results without KO to KEGG conversion
metacyc_daa_annotated_results_df <- pathway_annotation(pathway = "MetaCyc", daa_results_df = abundance_daa_results_df, ko_to_kegg = FALSE)

# Generate pathway heatmap
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset
feature_with_p_0.05 <- abundance_daa_results_df %>% filter(p_values < 0.05)
pd_antidep_heatmap <- pathway_heatmap(abundance = abundance_data_filtered %>% filter(pathway %in% feature_with_p_0.05$feature) %>% column_to_rownames("pathway"), metadata = PD_metadata, group = "Antidepressant_use")
pd_antidep_heatmap

## Generate Heatmap with Description ##
# Add description of the metabolic pathways to filtered abundance data, relocate the description to the first column
abundance_data_filtered_with_description <- cbind(abundance_data_filtered, description = metacyc_daa_annotated_results_df$description) %>%
  relocate(description)

# Remove the pathway column
abundance_data_filtered_with_description = abundance_data_filtered_with_description[,-2]

description_p_0.05 <- metacyc_daa_annotated_results_df %>% filter(p_values < 0.05)
pd_antidepressant_map <- pathway_heatmap(abundance = abundance_data_filtered_with_description %>% filter(description %in% description_p_0.05$description) %>% column_to_rownames("description"), 
                                         metadata = PD_metadata, 
                                         group = "Antidepressant_use")
pd_antidepressant_map
ggsave(filename = "pd_antidep_heatmap.png", pd_antidepressant_map, height = 11, width = 17)

# Generate pathway PCA plot
# Please change column_to_rownames() to the feature column if you are not using example dataset
# Please change group to "your_group_column" if you are not using example dataset
pd_antidepressant_pca <- pathway_pca(abundance = abundance_data_filtered %>% column_to_rownames("pathway"), metadata = PD_metadata, group = "Antidepressant_use")
pd_antidepressant_pca

ggsave(filename = "pd_antidep_pca.png", pd_antidepressant_pca, height = 8, width = 10)

res =  DEseq2_function(abundance_data_filtered,PD_metadata,"Antidepressant_use")
res$feature =rownames(res)
res_desc = inner_join(res,metacyc_daa_annotated_results_df, by = "feature")
res_desc = res_desc[, -c(8:13)]
View(res_desc)

sig_res = res_desc %>%
  filter(pvalue < 0.05)

sig_res <- sig_res[order(sig_res$log2FoldChange),]
pd_antidep_log <- ggplot(data = sig_res, aes(y = reorder(description, sort(as.numeric(log2FoldChange))), x= log2FoldChange, fill = pvalue))+
  geom_bar(stat = "identity")+
  theme_bw()+
  labs(x = "Log2FoldChange", y="Pathways")
pd_antidep_log

ggsave(filename = "pd_antidep_log.png", pd_antidep_log, height = 7, width = 9)
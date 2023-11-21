#### Load packages ####
library(phyloseq)
library(ape) 
library(tidyverse)
library(vegan)

#add column combining depression_binned and Disease status data
meta <- read.csv("parkinsons_export/parkinsons_metadata_new_edited.csv")
dep_PD_meta <- meta %>% unite("depression_binned_Disease", 
                              depression_binned, Disease, remove=F)

### Load data for phyloseq object ###
meta <- dep_PD_meta

otufp <- "parkinsons_export/feature-table.txt"
otu <- read_delim(file = otufp, delim = "\t", skip =1)

taxfp <- "parkinsons_export/taxonomy.tsv"
tax <- read_delim(taxfp, delim = "\t")

phylotreefp <- "parkinsons_export/tree.nwk"
phylotree <- read.tree(phylotreefp)

### Format OTU table ###
otu_mat <- as.matrix(otu[,-1])
rownames(otu_mat) <- otu$'#OTU ID'
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)

### Format sample metadata ###
samp_df <- as.data.frame(meta[,-1])
rownames(samp_df) <- meta$'X.SampleID'

#Make phyloseq sample data with sample_data() function
SAMP <- sample_data(samp_df)

### Formatting taxonomy ###
#convert taxon strings to a table with separate taxa rank columns 
tax_mat <- tax %>% select(-Confidence)%>%
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix() 

# Save everything except feature IDs 
tax_mat <- tax_mat[,-1]

# Make sampleids the rownames
rownames(tax_mat) <- tax$`Feature ID`

# Make taxa table
TAX <- tax_table(tax_mat)

### Make phyloseq object ###
parkinsons <- phyloseq(OTU, SAMP, TAX, phylotree)

### Analyze phyloseq object ###
#Remove non-bacterial sequences if any
parkinsons_filt <- subset_taxa(parkinsons,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")

#Remove ASVs with less than 5 counts total 
parkinsons_filtnolow <- filter_taxa(parkinsons_filt, function(x) sum(x)>5, prune = TRUE) 

# Remove samples with less than 100 reads
parkinsons_final <- prune_samples(sample_sums(parkinsons_filtnolow)>100, parkinsons_filtnolow)

# Rarefy samples for diversity metrics
rarecurve(t(as.data.frame(otu_table(parkinsons_final))), cex=0.1) 
parkinsons_rare <- rarefy_even_depth(parkinsons_final, rngseed = 1, sample.size = 3797) 

save(parkinsons_final, file = "parkinsons_dep_disease_final.RData")
save(parkinsons_rare, file = "parkinsons_dep_disease_rare.RData")

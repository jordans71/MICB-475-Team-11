#Import libraries
library(tidyverse)
library(phyloseq)
library(ape) 
library(vegan)
library(picante)
### Load all of the necessary files to make phyloseq object ###
#recreate phlyoseq object with new metadata containing binned depression and anxiety scores
metafp <- "parkinsons_metadata_new_edited.csv"
meta <- read_delim(metafp, delim=",")

otufp <- "feature-table.txt"
otu <- read_delim(file = otufp, delim="\t", skip=1)

taxfp <- "taxonomy.tsv"
tax <- read_delim(taxfp, delim="\t")

phylotreefp <- "tree.nwk"
phylotree <- read.tree(phylotreefp)

## Format OTU Table ##
otu_mat <- as.matrix(otu[,-1])

# Make first column (#OTU ID) the rownames of the new matrix
rownames(otu_mat) <- otu$`#OTU ID`

# Use the "otu_table" function to make an OTU table
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) 

## Format metadata ##
# Save everything except sampleid as new data frame
samp_df <- as.data.frame(meta[,-2])

# Make sampleids the rownames
rownames(samp_df)<- meta$'X.SampleID'

# Make phyloseq sample data with sample_data() function
SAMP <- sample_data(samp_df)

## Format taxonomy ##
# Convert taxon strings to a table with separate taxa rank columns
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

## Create phyloseq object ##
parkinsons <- phyloseq(OTU, SAMP, TAX, phylotree)

#check objects
otu_table(parkinsons)
sample_data(parkinsons)
tax_table(parkinsons)
phy_tree(parkinsons)


######### ANALYZE ##########
# Remove non-bacterial sequences, if any
parkinsons_filt <- subset_taxa(parkinsons,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")

# Remove samples with less than 100 reads
parkinsons_filt_nolow_samps <- prune_samples(sample_sums(parkinsons_filt)>100, parkinsons_filt)


### RAREFY###
rarecurve(t(as.data.frame(otu_table(parkinsons_filt_nolow_samps))), cex=0.1)
parkinsons_rare <- rarefy_even_depth(parkinsons_filt_nolow_samps, rngseed = 1, sample.size = 3797)

#to save rarefied phyloseq object, filtered phyloseq object 
save(parkinsons_filt_nolow_samps, file="parkinsons2_filt_nolow_samps.RData")
save(parkinsons_rare, file="parkinsons_edited_rare.RData")

View(sample_data(parkinsons_rare))

#create specific phyloseq object for each aim
# Remove samples where anxiety is na
parkinsons_final_anxiety <- subset_samples(parkinsons_rare, !is.na(anxiety_binned))
View(sample_data(parkinsons_final_anxiety))
save(parkinsons_final_anxiety, file="parkinsons_final_anxiety.RData")

#remove samples where sleep is na
parkinsons_final_sleep <- subset_samples(parkinsons_rare, !is.na(Sleep_problems))
View(sample_data(parkinsons_final_sleep))
save(parkinsons_final_sleep, file="parkinsons_final_sleep.RData")

#remove samples where depression is na
parkinsons_final_depression <- subset_samples(parkinsons_rare, !is.na(depression_binned))
View(sample_data(parkinsons_final_depression))
save(parkinsons_final_depression, file="parkinsons_final_depression.RData")

############ Anxiety ###########
### Load anxiety phyloseq Object 
load("1_make_phyloseq_objects/parkinsons_final_anxiety.RData")

#Filter out PD patients 
PD_patients <- subset_samples(parkinsons_final_anxiety, `Disease` == "PD")

#Filter out control
Ctrl_patients <- subset_samples(parkinsons_final_anxiety, `Disease` == "Control")

gg_richness_PD_anxiety <- plot_richness(PD_patients, x = "anxiety_binned", measures = c("shannon")) + 
  xlab("PD Anxiety") + geom_boxplot()
gg_richness_Ctrl_anxiety <- plot_richness(Ctrl_patients, x = "anxiety_binned", measures = c("shannon")) + 
  xlab("Control Anxiety") + geom_boxplot()

ggsave("PD_anxiety_shannon.png"
       , gg_richness_PD_anxiety 
       , height=4, width=6)
ggsave("ctrl_anxiety_shannon.png"
       , gg_richness_Ctrl_anxiety
       , height=4, width=6)
anxiety_shannon <- plot_grid(gg_richness_PD_anxiety, gg_richness_Ctrl_anxiety, labels = c('C', 'D'))
anxiety_shannon

### Load depression phyloseq Object 
load("1_make_phyloseq_objects/parkinsons_final_depression.RData")

#Filter out PD patients 
PD_patients <- subset_samples(parkinsons_final_depression, `Disease` == "PD")

#Filter out control
Ctrl_patients <- subset_samples(parkinsons_final_depression, `Disease` == "Control")

gg_richness_PD_depression <- plot_richness(PD_patients, x = "depression_binned", measures = c("shannon")) + 
  xlab("PD Depression") + geom_boxplot()
gg_richness_Ctrl_depression <- plot_richness(Ctrl_patients, x = "depression_binned", measures = c("shannon")) + 
  xlab("Control Depression") + geom_boxplot()

ggsave("PD_depression_shannon.png"
       , gg_richness_PD_depression 
       , height=4, width=6)
ggsave("ctrl_depression_shannon.png"
       , gg_richness_Ctrl_depression
       , height=4, width=6)
depression_shannon <- plot_grid(gg_richness_PD_depression, gg_richness_Ctrl_depression, labels = c('A', 'B'))
depression_shannon

### Load depression phyloseq Object 
load("1_make_phyloseq_objects/parkinsons_final_sleep_problems.RData")

#Filter out PD patients 
PD_patients <- subset_samples(parkinsons_final_sleep_problems, `Disease` == "PD")

#Filter out control
Ctrl_patients <- subset_samples(parkinsons_final_sleep_problems, `Disease` == "Control")

gg_richness_PD_sleep <- plot_richness(PD_patients, x = "Sleep_problems", measures = c("shannon")) + 
  xlab("PD Sleep Problems") + geom_boxplot()
gg_richness_Ctrl_sleep <- plot_richness(Ctrl_patients, x = "Sleep_problems", measures = c("shannon")) + 
  xlab("Control Sleep Problems") + geom_boxplot()

ggsave("PD_sleep_shannon.png"
       , gg_richness_PD_sleep
       , height=4, width=6)
ggsave("ctrl_sleep_shannon.png"
       , gg_richness_Ctrl_sleep
       , height=4, width=6)
sleep_shannon <- plot_grid(gg_richness_PD_sleep, gg_richness_Ctrl_sleep, labels = c('E', 'F'))
sleep_shannon
##### PD #####

### Load libraries ###
library(phyloseq)
library(ape) 
library(tidyverse)
library(vegan)
### Load depression phyloseq Object 
### Load all of the necessary files to make phyloseq object ###
#recreate phlyoseq object with new metadata containing binned depression and anxiety scores
metafp <- "5_aim4_disease/parkinsons_metadata_new_edited.csv"
meta <- read_delim(metafp, delim=",")

otufp <- "5_aim4_disease/feature-table.txt"
otu <- read_delim(file = otufp, delim="\t", skip=1)

taxfp <- "5_aim4_disease/taxonomy.tsv"
tax <- read_delim(taxfp, delim="\t")

phylotreefp <- "5_aim4_disease/tree.nwk"
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

gg_richness_PD <- plot_richness(parkinsons_rare, x = "Disease", measures = c("shannon")) + 
  xlab("Disease Status") + geom_boxplot()

ggsave("PD_shannon.png"
       , gg_richness_PD
       , height=4, width=6)
PD_shannon <- plot_grid(gg_richness_PD,  labels = c('G'))
PD_shannon

dep_anxiety_sleep_disease_together_shannon <- grid.arrange(depression_shannon, anxiety_shannon,sleep_shannon, PD_shannon, ncol = 1)
dep_anxiety_sleep_disease_together_shannon
ggsave("Shannon.png"
       , dep_anxiety_sleep_disease_together_shannon
       , height=15, width=10)



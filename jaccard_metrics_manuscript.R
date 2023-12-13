#load libraries
library(phyloseq)
library(ape) 
library(tidyverse)
library(vegan)
library(cowplot)

############ Anxiety ###########
### Load anxiety phyloseq Object 
load("1_make_phyloseq_objects/parkinsons_final_anxiety.RData")

#Filter out PD patients 
PD_patients <- subset_samples(parkinsons_final_anxiety, `Disease` == "PD")

#Filter out control
Ctrl_patients <- subset_samples(parkinsons_final_anxiety, `Disease` == "Control")

###Data frame 
#PD
samp_dat_wdiv_PD <- data.frame(sample_data(PD_patients), estimate_richness(PD_patients))

#Control
samp_dat_wdiv_Ctrl <- data.frame(sample_data(Ctrl_patients), estimate_richness(Ctrl_patients))

## Jaccard ## 
#PD patients
jac_dm <- distance(PD_patients, method = "jaccard", binary = T)
pcoa_jac_PD <- ordinate(PD_patients, method = "NMDS", distance = jac_dm)
PD_anxiety_jac <- plot_ordination(PD_patients, pcoa_jac_PD, color = "anxiety_binned") +
  labs(col = "") + theme_bw() + stat_ellipse(level = 0.95) +
  ggtitle("Jaccard PD Anxiety") + theme(plot.title = element_text(hjust = 0.5)) 
PD_anxiety_jac
ggsave("PD_anxiety_jac_pcoa.png"
       , PD_anxiety_jac
       , height=4, width=6)

adonis2(jac_dm ~ `anxiety_binned`, data = samp_dat_wdiv_PD)

#Healthy Controls
jac_dm_ctrl <- distance(Ctrl_patients, method = "jaccard", binary = T)
pcoa_jac_ctrl <- ordinate(Ctrl_patients, method = "NMDS", distance = jac_dm_ctrl)
ctrl_anxiety_jac <- plot_ordination(Ctrl_patients, pcoa_jac_ctrl, color = "anxiety_binned") + 
  labs(col = "Anxiety") + theme_bw() + stat_ellipse(level = 0.95) +
  ggtitle("Jaccard Control Anxiety") + theme(plot.title = element_text(hjust = 0.5)) 
ctrl_anxiety_jac

ggsave("ctrl_anxiety_jac_pcoa.png"
       , ctrl_anxiety_jac
       , height=4, width=6)
adonis2(jac_dm_ctrl ~ `anxiety_binned`, data = samp_dat_wdiv_Ctrl)

anxiety_jac <- plot_grid(PD_anxiety_jac, ctrl_anxiety_jac, labels = c('C', 'D'))
anxiety_jac

################ Depression ################
### Load depression phyloseq Object 
load("1_make_phyloseq_objects/parkinsons_final_depression.RData")

#Filter out PD patients 
PD_patients <- subset_samples(parkinsons_final_depression, `Disease` == "PD")

#Filter out control
Ctrl_patients <- subset_samples(parkinsons_final_depression, `Disease` == "Control")

###Data frame 
#PD
samp_dat_wdiv_PD <- data.frame(sample_data(PD_patients), estimate_richness(PD_patients))

#Control
samp_dat_wdiv_Ctrl <- data.frame(sample_data(Ctrl_patients), estimate_richness(Ctrl_patients))

## Jaccard ## 
#PD patients
jac_dm <- distance(PD_patients, method = "jaccard", binary = T)
pcoa_jac_PD <- ordinate(PD_patients, method = "NMDS", distance = jac_dm)
PD_depression_jac <- plot_ordination(PD_patients, pcoa_jac_PD, color = "depression_binned") +
  labs(col = "") + theme_bw() + stat_ellipse(level = 0.95) +
  ggtitle("Jaccard PD Depression") + theme(plot.title = element_text(hjust = 0.5)) +
  guides(fill = "none") 
PD_depression_jac

ggsave("PD_depression_jac_pcoa.png"
       , PD_depression_jac
       , height=4, width=6)

adonis2(jac_dm ~ `depression_binned`, data = samp_dat_wdiv_PD)

#Healthy Controls
jac_dm_ctrl <- distance(Ctrl_patients, method = "jaccard", binary = T)
pcoa_jac_ctrl <- ordinate(Ctrl_patients, method = "NMDS", distance = jac_dm_ctrl)
ctrl_depression_jac <- plot_ordination(Ctrl_patients, pcoa_jac_ctrl, color = "depression_binned") + 
  labs(col = "Depression") + theme_bw() + stat_ellipse(level = 0.95) +
  ggtitle("Jaccard Control Depression") + theme(plot.title = element_text(hjust = 0.5)) 
ctrl_depression_jac
ggsave("ctrl_depression_jac_pcoa.png"
       , ctrl_depression_jac
       , height=4, width=6)
adonis2(jac_dm_ctrl ~ `depression_binned`, data = samp_dat_wdiv_Ctrl)

depression_jac <- plot_grid(PD_depression_jac, ctrl_depression_jac, labels = c('A', 'B'))
depression_jac


####Sleep problems ####
load("1_make_phyloseq_objects/parkinsons_final_sleep_problems.RData")

#Filter out PD patients 
PD_patients <- subset_samples(parkinsons_final_sleep_problems , `Disease` == "PD")

#Filter out control
Ctrl_patients <- subset_samples(parkinsons_final_sleep_problems , `Disease` == "Control")

###Data frame 
#PD
samp_dat_wdiv_PD <- data.frame(sample_data(PD_patients), estimate_richness(PD_patients))

#Control
samp_dat_wdiv_Ctrl <- data.frame(sample_data(Ctrl_patients), estimate_richness(Ctrl_patients))

## Jaccard ## 
#PD patients
jac_dm <- distance(PD_patients, method = "jaccard", binary = T)
pcoa_jac_PD <- ordinate(PD_patients, method = "NMDS", distance = jac_dm)
PD_Sleep_problems_jac <- plot_ordination(PD_patients, pcoa_jac_PD, color = "Sleep_problems") +
  labs(col = "") + theme_bw() + stat_ellipse(level = 0.95) +
  ggtitle("Jaccard PD Sleep") + theme(plot.title = element_text(hjust = 0.5)) +
  guides(fill = "none") 
PD_Sleep_problems_jac

ggsave("PD_Sleep_problems_jac_pcoa.png"
       , PD_Sleep_problems_jac
       , height=4, width=6)

adonis2(jac_dm ~ `Sleep_problems`, data = samp_dat_wdiv_PD)

#Healthy Controls
jac_dm_ctrl <- distance(Ctrl_patients, method = "jaccard", binary = T)
pcoa_jac_ctrl <- ordinate(Ctrl_patients, method = "NMDS", distance = jac_dm_ctrl)
ctrl_Sleep_problems_jac <- plot_ordination(Ctrl_patients, pcoa_jac_ctrl, color = "Sleep_problems") + 
  labs(col = "Sleep problems") + theme_bw() + stat_ellipse(level = 0.95) +
  ggtitle("Jaccard Control Sleep") + theme(plot.title = element_text(hjust = 0.5)) 
ctrl_Sleep_problems_jac

ggsave("ctrl_Sleep_problems _jac_pcoa.png"
       , ctrl_depression_jac
       , height=4, width=6)
adonis2(jac_dm_ctrl ~ `Sleep_problems`, data = samp_dat_wdiv_Ctrl)

Sleep_problems_jac <- plot_grid(PD_Sleep_problems_jac, ctrl_Sleep_problems_jac, labels = c('E', 'F'))
Sleep_problems_jac
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


###Data frame 
samp_dat_wdiv <- data.frame(sample_data(parkinsons_rare), estimate_richness(parkinsons_rare))

## Jaccard ## 
#PD patients
jac_dm <- distance(parkinsons_rare, method = "jaccard", binary = T)
pcoa_jac_PD <- ordinate(parkinsons_rare, method = "NMDS", distance = jac_dm)
PD_jac <- plot_ordination(parkinsons_rare, pcoa_jac_PD, color = "Disease") +
  labs(col = "Disease Status") + theme_bw() + stat_ellipse(level = 0.95) +
  ggtitle("Jaccard Disease") + theme(plot.title = element_text(hjust = 0.5)) +
  guides(fill = "none") 
PD_jac <- plot_grid(PD_jac, labels = c('G'))

ggsave("PD_jac_pcoa.png"
       , PD_jac
       , height=4, width=6)
PD_jac




library(gridExtra)
dep_anxiety_sleep_disease_together <- grid.arrange(depression_jac, anxiety_jac,Sleep_problems_jac, PD_jac, ncol = 1)
dep_anxiety_sleep_disease_together
ggsave("Jac_Pcoa.png"
       , dep_anxiety_sleep_disease_together
       , height=15, width=10)

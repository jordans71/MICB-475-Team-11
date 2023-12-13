### Load libraries ###
library(phyloseq)
library(ape) 
library(tidyverse)
library(vegan)
### Load depression phyloseq Object 
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



### Beta Diversity ###
## Jaccard ## 
jac_dm <- distance(parkinsons_rare, method = "jaccard", binary = TRUE)
pcoa_jac <- ordinate(parkinsons_rare, method = "NMDS", distance = jac_dm)
gg_jac_pcoa <- plot_ordination(parkinsons_rare, pcoa_jac, color = "Disease") +
  labs(col = "Disease Status") + theme_bw() + stat_ellipse(level = 0.95) +
  ggtitle("Jaccard_Disease") + theme(plot.title = element_text(hjust = 0.5)) +
  ggside::geom_xsideboxplot(aes(fill = Disease, y = Disease), 
                            orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = Disease, x = Disease), 
                            orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void()
  
gg_jac_pcoa

ggsave("Disease_jac_pcoa.png"
       , gg_jac_pcoa
       , height=4, width=6)

adonis2(jac_dm ~ Disease, data = samp_dat_wdiv)

## bray curtis ##
bray_dm <- distance(parkinsons_rare, method = "bray")
pcoa_bray <- ordinate(parkinsons_rare, method = "PCoA", distance = bray_dm)
gg_bray_pcoa <- plot_ordination(parkinsons_rare, pcoa_bray, color = "Disease") +
  labs(col = "Disease Status") + theme_bw() + stat_ellipse(level = 0.95) +
  ggtitle("Bray_Curtis_Disease") + theme(plot.title = element_text(hjust = 0.5)) +
  ggside::geom_xsideboxplot(aes(fill = Disease, y = Disease), 
                            orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = Disease, x = Disease), 
                            orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void() 
gg_bray_pcoa

ggsave("Disease_bray_pcoa.png"
       , gg_bray_pcoa
       , height=4, width=6)

adonis2(bray_dm ~ Disease, data = samp_dat_wdiv)

## unweighted unifrac ##
unifrac_dm <- distance(parkinsons_rare, method = "unifrac")
pcoa_unifrac <- ordinate(parkinsons_rare, method = "PCoA", distance = unifrac_dm)
gg_unifrac_pcoa <- plot_ordination(parkinsons_rare, pcoa_unifrac, color = "Disease") +
  labs(col = "Disease Status") + theme_bw() + stat_ellipse(level = 0.95) +
  ggtitle("Unweighted_Unifrac_Disease") + theme(plot.title = element_text(hjust = 0.5))+
  ggside::geom_xsideboxplot(aes(fill = Disease, y = Disease), 
                            orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = Disease, x = Disease), 
                            orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void() 
gg_unifrac_pcoa

ggsave("Disease_unifrac_pcoa.png"
       , gg_unifrac_pcoa
       , height=4, width=6)

adonis2(unifrac_dm ~ Disease, data = samp_dat_wdiv)

## weighted_unifrac ##
w_unifrac_dm <- distance(parkinsons_rare, method ="wunifrac")
pcoa_w_unifrac <- ordinate(parkinsons_rare, method="PCoA", distance=w_unifrac_dm)
gg_wunifrac_pcoa <- plot_ordination(parkinsons_rare, pcoa_w_unifrac, color = "Disease") +
  labs(col = "Disease Status") + theme_bw() + stat_ellipse(level = 0.95) +
  ggtitle("Weighted_Unifrac_Disease") + theme(plot.title = element_text(hjust = 0.5)) +
  ggside::geom_xsideboxplot(aes(fill = Disease, y = Disease), 
                            orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = Disease, x = Disease), 
                            orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void() 
gg_wunifrac_pcoa

ggsave("Disease_wunifrac_pcoa.png"
       , gg_wunifrac_pcoa
       , height=4, width=6)

adonis2(w_unifrac_dm ~ Disease, data = samp_dat_wdiv)

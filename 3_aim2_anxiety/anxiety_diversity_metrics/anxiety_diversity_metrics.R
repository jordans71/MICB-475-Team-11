### Load libraries ###
library(phyloseq)
library(ape) 
library(tidyverse)
library(vegan)

#add column combining anxiety_binned and Disease status data
meta <- read.csv("parkinsons_metadata_new_edited.csv")

anx_Disease_meta <- meta %>% unite("anxiety_binned_Disease", 
                                   anxiety_binned, Disease, remove=F)
### Load data for phyloseq object ###
meta <- anx_Disease_meta

otufp <- "feature-table.txt"
otu <- read_delim(file = otufp, delim = "\t", skip =1)

taxfp <- "taxonomy.tsv"
tax <- read_delim(taxfp, delim = "\t")

phylotreefp <- "tree.nwk"
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

save(parkinsons_final, file = "parkinsons_anxiety_disease_final.RData")
save(parkinsons_rare, file = "parkinsons_anxiety_disease_rare.RData")

# Remove samples where anxiety is na
parkinsons_final_anxiety <- subset_samples(parkinsons_rare, !is.na(anxiety_binned))
View(sample_data(parkinsons_final_anxiety))
save(parkinsons_final_anxiety, file="parkinsons_final_anxiety.RData")

load("parkinsons_final_anxiety.RData")

### Alpha Diversity ###
gg_richness <- plot_richness(parkinsons_final_anxiety, x = "anxiety_binned_Disease") + 
  xlab("anxiety_PD_Status") + geom_boxplot()
gg_richness

###Statistical analysis 
samp_dat_wdiv <- data.frame(sample_data(parkinsons_final_anxiety), estimate_richness(parkinsons_final_anxiety))

ggplot(samp_dat_wdiv) + geom_boxplot(aes(x=anxiety_binned, y=Shannon)) +
  facet_grid(~factor(`Disease`, levels=c("PD","Control")))

# run the 2-way ANOVA on Shannon diversity 
ml_anxiety_sub <- lm(Shannon ~ `anxiety_binned`*`Disease`, data=samp_dat_wdiv)
summary(aov(ml_anxiety_sub))
TukeyHSD(aov(ml_anxiety_sub))

### Linear models ####
# Linear models are identical to ANOVA when predictors are categorical
ggplot(samp_dat_wdiv) + geom_boxplot(aes(x=anxiety_binned, y=Shannon)) +
  facet_grid(~factor(`Disease`))
summary(ml_anxiety_sub)


### Beta Diversity ###
## Jaccard ## 
#PD patients
jac_dm <- distance(parkinsons_final_anxiety, method = "jaccard", binary = T)
pcoa_jac_PD <- ordinate(parkinsons_final_anxiety, method = "NMDS", distance = jac_dm)
PD_anxiety_jac <- plot_ordination(parkinsons_final_anxiety, pcoa_jac_PD, color = "anxiety_binned") +
  labs(col = "Anxiety Level") + theme_bw() + stat_ellipse(level = 0.95) +
  ggside::geom_xsideboxplot(aes(fill = anxiety_binned, y = anxiety_binned), orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = anxiety_binned, x = anxiety_binned), orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void() 
PD_anxiety_jac

#Healthy Controls
jac_dm_ctrl <- distance(parkinsons_final_anxiety, method = "jaccard", binary = T)
pcoa_jac_ctrl <- ordinate(parkinsons_final_anxiety, method = "NMDS", distance = jac_dm_ctrl)
ctrl_anxiety_jac <- plot_ordination(parkinsons_final_anxiety, pcoa_jac_ctrl, color = "anxiety_binned") + 
  labs(col = "Anxiety Level") + theme_bw() + stat_ellipse(level = 0.95) +
  ggside::geom_xsideboxplot(aes(fill = anxiety_binned, y = anxiety_binned), orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = anxiety_binned, x = anxiety_binned), orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void() 
ctrl_anxiety_jac

##Bray-Curtis ##
# PD patients
bc_dm <- distance(parkinsons_final_anxiety, method="bray")
pcoa_bc_PD <- ordinate(parkinsons_final_anxiety, method="PCoA", distance=bc_dm)
PD_anxiety_bray <- plot_ordination(parkinsons_final_anxiety, pcoa_bc_PD, color = "anxiety_binned") + 
  labs(col = "anxiety level") + theme_bw() + stat_ellipse(level = 0.95) +
  ggside::geom_xsideboxplot(aes(fill = anxiety_binned, y = anxiety_binned), orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = anxiety_binned, x = anxiety_binned), orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void()
PD_anxiety_bray

## Healthy controls
bc_dm_ctrl <- distance(parkinsons_final_anxiety, method="bray")
pcoa_bc_ctrl <- ordinate(parkinsons_final_anxiety, method="PCoA", distance=bc_dm_ctrl)
ctrl_anxiety_bray <- plot_ordination(parkinsons_final_anxiety, pcoa_bc_ctrl, color = "anxiety_binned") + 
  labs(col = "anxiety level") + theme_bw() + stat_ellipse(level = 0.95) +
  ggside::geom_xsideboxplot(aes(fill = anxiety_binned, y = anxiety_binned), orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = anxiety_binned, x = anxiety_binned), orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void()
ctrl_anxiety_bray

## Unweighted Unifrac ##
#PD patients 
unifrac_dm_PD <- distance(parkinsons_final_anxiety, method = "unifrac")
pcoa_unifrac_PD <- ordinate(parkinsons_final_anxiety, method = "PCoA", distance = unifrac_dm_PD)
PD_gg_unifrac_pcoa <- plot_ordination(parkinsons_final_anxiety, pcoa_unifrac_PD, color = "anxiety_binned") +
  labs(col = "Anxiety Level") + theme_bw() + stat_ellipse(level = 0.95) +
  ggtitle("Unweighted Unifrac") + theme(plot.title = element_text(hjust = 0.5)) +
  ggside::geom_xsideboxplot(aes(fill = anxiety_binned, y = anxiety_binned), orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = anxiety_binned, x = anxiety_binned), orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void() 
PD_gg_unifrac_pcoa

#Healthy Controls
unifrac_dm_ctrl <- distance(parkinsons_final_anxiety, method = "unifrac")
pcoa_unifrac_ctrl <- ordinate(parkinsons_final_anxiety, method = "PCoA", distance = unifrac_dm_ctrl)
ctrl_gg_unifrac_pcoa <- plot_ordination(parkinsons_final_anxiety, pcoa_unifrac_ctrl, color = "anxiety_binned") +
  labs(col = "Anxiety Level") + theme_bw() + stat_ellipse(level = 0.95) +
  ggtitle("Unweighted Unifrac") + theme(plot.title = element_text(hjust = 0.5)) +
  ggside::geom_xsideboxplot(aes(fill = anxiety_binned, y = anxiety_binned), orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = anxiety_binned, x = anxiety_binned), orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void() 
ctrl_gg_unifrac_pcoa

## Weighted Unifrac ##
#PD patients 
w_unifrac_dm_PD <- distance(parkinsons_final_anxiety, method ="wunifrac")
pcoa_w_unifrac_PD <- ordinate(parkinsons_final_anxiety, method="PCoA", distance=w_unifrac_dm_PD)
gg_wunifrac_pcoa_PD <- plot_ordination(parkinsons_final_anxiety, pcoa_w_unifrac_PD, color = "anxiety_binned") +
  labs(col = "Anxiety Status") + theme_bw() + stat_ellipse(level = 0.95) +
  ggtitle("Weighted Unifrac") + theme(plot.title = element_text(hjust = 0.5)) +
  ggside::geom_xsideboxplot(aes(fill = anxiety_binned, y = anxiety_binned), orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = anxiety_binned, x = anxiety_binned), orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void() 
gg_wunifrac_pcoa_PD

#Healthy Controls
w_unifrac_dm_ctrl <- distance(parkinsons_final_anxiety, method ="wunifrac")
pcoa_w_unifrac_ctrl <- ordinate(parkinsons_final_anxiety, method="PCoA", distance=w_unifrac_dm_ctrl)
gg_wunifrac_pcoa_ctrl <- plot_ordination(parkinsons_final_anxiety, pcoa_w_unifrac_ctrl, color = "anxiety_binned") +
  labs(col = "Anxiety Status") + theme_bw() + stat_ellipse(level = 0.95) +
  ggtitle("Weighted Unifrac") + theme(plot.title = element_text(hjust = 0.5)) +
  ggside::geom_xsideboxplot(aes(fill = anxiety_binned, y = anxiety_binned), orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = anxiety_binned, x = anxiety_binned), orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void() 
gg_wunifrac_pcoa_ctrl

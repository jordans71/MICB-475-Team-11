### Load libraries ###
library(phyloseq)
library(ape) 
library(tidyverse)
library(vegan)

#add column combining sleep_problems and Disease status data
meta <- read.csv("parkinsons_metadata_new_edited.csv")

sleep_problems_Disease_meta <- meta %>% unite("Sleep_problems_binned_Disease", 
                                              Sleep_problems, Disease, remove=F)
### Load data for phyloseq object ###
meta <- sleep_problems_Disease_meta

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

save(parkinsons_final, file = "parkinsons_sleep_problems_disease_final.RData")
save(parkinsons_rare, file = "parkinsons_sleep_problems_disease_rare.RData")

# Remove samples where anxiety is na
parkinsons_final_sleep_problems <- subset_samples(parkinsons_rare, Sleep_problems != "")
View(sample_data(parkinsons_final_sleep_problems))
save(parkinsons_final_sleep_problems, file="parkinsons_final_sleep_problems.RData")

load("parkinsons_final_sleep_problems.RData")

### Alpha Diversity ###
gg_richness <- plot_richness(parkinsons_final_sleep_problems, x = "Sleep_problems_binned_Disease") + 
  xlab("sleep_PD_Status") + geom_boxplot()
gg_richness

###Statistical analysis 
samp_dat_wdiv <- data.frame(sample_data(parkinsons_final_sleep_problems), estimate_richness(parkinsons_final_sleep_problems))

ggplot(samp_dat_wdiv) + geom_boxplot(aes(x=Sleep_problems, y=Shannon)) +
  facet_grid(~factor(`Disease`, levels=c("PD","Control")))

# run the 2-way ANOVA on Shannon diversity 
ml_sleep_sub <- lm(Shannon ~ `Sleep_problems`*`Disease`, data=samp_dat_wdiv)
summary(aov(ml_sleep_sub))
TukeyHSD(aov(ml_sleep_sub))

### Linear models ####
# Linear models are identical to ANOVA when predictors are categorical
ggplot(samp_dat_wdiv) + geom_boxplot(aes(x=Sleep_problems, y=Shannon)) +
  facet_grid(~factor(`Disease`))
summary(ml_sleep_sub)

### Beta Diversity ###
## Jaccard ## 
#PD patients
jac_dm <- distance(PD_patients, method = "jaccard", binary = T)
pcoa_jac_PD <- ordinate(PD_patients, method = "NMDS", distance = jac_dm)
PD_sleep_jac <- plot_ordination(PD_patients, pcoa_jac_PD, color = "Sleep_problems") +
  labs(col = "Sleep Problems") + theme_bw() + stat_ellipse(level = 0.95) +
  ggside::geom_xsideboxplot(aes(fill = Sleep_problems, y = Sleep_problems), orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = Sleep_problems, x = Sleep_problems), orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void() +
  ggtitle("PD Sleep Jaccard") + theme(plot.title = element_text(hjust = 0.5))
PD_sleep_jac

#Healthy Controls
jac_dm_ctrl <- distance(ctrl_patients, method = "jaccard", binary = T)
pcoa_jac_ctrl <- ordinate(ctrl_patients, method = "NMDS", distance = jac_dm_ctrl)
ctrl_sleep_jac <- plot_ordination(ctrl_patients, pcoa_jac_ctrl, color = "Sleep_problems") + 
  labs(col = "Sleep_problems") + theme_bw() + stat_ellipse(level = 0.95) +
  ggside::geom_xsideboxplot(aes(fill = Sleep_problems, y = Sleep_problems), orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = Sleep_problems, x = Sleep_problems), orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void() +
  ggtitle("Healthy Sleep Jaccard") + theme(plot.title=element_text(hjust = 0.5))
ctrl_sleep_jac

##Bray-Curtis ##
# PD patients
bc_dm <- distance(PD_patients, method="bray")
pcoa_bc <- ordinate(PD_patients, method="PCoA", distance=bc_dm)
PD_sleep_bray <- plot_ordination(PD_patients, pcoa_bc, color = "Sleep_problems") + 
  labs(col = "Sleep_problems") + theme_bw() + stat_ellipse(level = 0.95) +
  ggside::geom_xsideboxplot(aes(fill = Sleep_problems, y = Sleep_problems), orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = Sleep_problems, x = Sleep_problems), orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void() +
  ggtitle("PD Sleep Bray") + theme(plot.title = element_text(hjust=0.5))
PD_sleep_bray

## Healthy controls
bc_dm_ctrl <- distance(ctrl_patients, method="bray")
pcoa_bc_ctrl <- ordinate(ctrl_patients, method="PCoA", distance=bc_dm_ctrl)
ctrl_sleep_bray <- plot_ordination(ctrl_patients, pcoa_bc_ctrl, color = "Sleep_problems") + 
  labs(col = "Sleep_problems") + theme_bw() + stat_ellipse(level = 0.95) +
  ggside::geom_xsideboxplot(aes(fill = Sleep_problems, y = Sleep_problems), orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = Sleep_problems, x = Sleep_problems), orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void() +
  ggtitle("Healthy Sleep Bray") + theme(plot.title = element_text(hjust=0.5))
ctrl_sleep_bray

## Unweighted Unifrac ##
#PD patients 
unifrac_dm_PD <- distance(PD_patients, method = "unifrac")
pcoa_unifrac_PD <- ordinate(PD_patients, method = "PCoA", distance = unifrac_dm_PD)
gg_unifrac_pcoa_PD <- plot_ordination(PD_patients, pcoa_unifrac_PD, color = "Sleep_problems") +
  labs(col = "Sleep_problems") + theme_bw() + stat_ellipse(level = 0.95) +
  ggtitle("PD Sleep Unweighted Unifrac") + theme(plot.title = element_text(hjust = 0.5)) +
  ggside::geom_xsideboxplot(aes(fill = Sleep_problems, y = Sleep_problems), orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = Sleep_problems, x = Sleep_problems), orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void() 
gg_unifrac_pcoa_PD

#Healthy Controls
unifrac_dm_ctrl <- distance(ctrl_patients, method = "unifrac")
pcoa_unifrac_ctrl <- ordinate(ctrl_patients, method = "PCoA", distance = unifrac_dm_ctrl)
gg_unifrac_pcoa_ctrl <- plot_ordination(ctrl_patients, pcoa_unifrac_ctrl, color = "Sleep_problems") +
  labs(col = "Sleep_problems") + theme_bw() + stat_ellipse(level = 0.95) +
  ggtitle("Healthy Sleep Unweighted Unifrac") + theme(plot.title = element_text(hjust = 0.5)) +
  ggside::geom_xsideboxplot(aes(fill = Sleep_problems, y = Sleep_problems), orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = Sleep_problems, x = Sleep_problems), orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void() 
gg_unifrac_pcoa_ctrl

## Weighted Unifrac ##
#PD patients 
w_unifrac_dm_PD <- distance(PD_patients, method ="wunifrac")
pcoa_w_unifrac_PD <- ordinate(PD_patients, method="PCoA", distance=w_unifrac_dm_PD)
gg_wunifrac_pcoa_PD <- plot_ordination(PD_patients, pcoa_w_unifrac_PD, color = "Sleep_problems") +
  labs(col = "Sleep_problems") + theme_bw() + stat_ellipse(level = 0.95) +
  ggtitle("PD Sleep Weighted Unifrac") + theme(plot.title = element_text(hjust = 0.5)) +
  ggside::geom_xsideboxplot(aes(fill = Sleep_problems, y = Sleep_problems), orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = Sleep_problems, x = Sleep_problems), orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void() 
gg_wunifrac_pcoa_PD

#Healthy Controls
w_unifrac_dm_ctrl <- distance(ctrl_patients, method ="wunifrac")
pcoa_w_unifrac_ctrl <- ordinate(ctrl_patients, method="PCoA", distance=w_unifrac_dm_ctrl)
gg_wunifrac_pcoa_ctrl <- plot_ordination(ctrl_patients, pcoa_w_unifrac_ctrl, color = "Sleep_problems") +
  labs(col = "Sleep_problems") + theme_bw() + stat_ellipse(level = 0.95) +
  ggtitle("Healthy Sleep Weighted Unifrac") + theme(plot.title = element_text(hjust = 0.5)) +
  ggside::geom_xsideboxplot(aes(fill = Sleep_problems, y = Sleep_problems), orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = Sleep_problems, x = Sleep_problems), orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void() 
gg_wunifrac_pcoa_ctrl

sleep_jac <- grid.arrange(PD_anxiety_jac, PD_sleep_bray,ctrl_sleep_jac, ctrl_sleep_bray, ncol =2)
ggsave(filename = "sleep_bray_jac.png", sleep_jac)

sleep_unifrac <- grid.arrange(gg_unifrac_pcoa_PD, gg_wunifrac_pcoa_PD,gg_unifrac_pcoa_ctrl, gg_wunifrac_pcoa_ctrl, ncol = 2)
ggsave(filename = "sleep_unifrac.png", sleep_unifrac)

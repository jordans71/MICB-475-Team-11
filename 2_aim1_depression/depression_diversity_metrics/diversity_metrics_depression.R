library(phyloseq)
library(ape)
library(tidyverse)
library(picante)
library(gridExtra)
library(ggplot2)
library(ggside)

### Load in RData ###
load("parkinsons_final_depression.RData")

### Alpha Diversity ###
## Depression Status Only ##
gg_richness <- plot_richness(parkinsons_final_depression, x = "depression_binned") + 
  xlab("Depression_Status") + geom_boxplot()
gg_richness

ggsave(filename = "plot_richness.png", gg_richness, height = 4, width = 10)

## Depression and PD Status ##
gg_richness_dep_PD <- plot_richness(parkinsons_final_depression, x = "depression_binned_Disease") + 
  xlab("Depression_PD_Status") + geom_boxplot()
gg_richness_dep_PD

ggsave(filename = "plot_richness_depPD.png", gg_richness_dep_PD, height = 4, width = 10)

### Phylogenetic diversity ###
#calculate Faith's phylogenetic diversity as PD
phylo_dist <- pd(t(otu_table(parkinsons_final_depression)), phy_tree(parkinsons_final_depression), 
                 include.root = F)
sample_data(parkinsons_final_depression)$PD <- phylo_dist$PD
plot.pd <- ggplot(sample_data(parkinsons_final_depression), aes(depression_binned, PD)) + geom_boxplot() + 
  xlab("Subject #") + ylab("Phylogenetic Diversity")
plot.pd

### Beta Diversity ###
## Jaccard ##
jac_dm <- distance(parkinsons_final_depression, method = "jaccard", binary = TRUE)
pcoa_jac <- ordinate(parkinsons_final_depression, method = "NMDS", distance = jac_dm)
gg_jac_pcoa <- plot_ordination(parkinsons_final_depression, pcoa_jac, color = "depression_binned_Disease") +
  labs(col = "Depression and Disease Status") + theme_bw() + stat_ellipse(level = 0.95) +
  ggside::geom_xsideboxplot(aes(fill = depression_binned_Disease, y = depression_binned_Disease), 
                            orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = depression_binned_Disease, x = depression_binned_Disease), 
                            orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void() +
  ggtitle("Jaccard") + theme(plot.title = element_text(hjust=0.5))
gg_jac_pcoa

ggsave("jaccard_pcoa.png"
       , gg_jac_pcoa
       , height=4, width=5)

## bray curtis ##
bray_dm <- distance(parkinsons_final_depression, method = "bray")
pcoa_bray <- ordinate(parkinsons_final_depression, method = "PCoA", distance = bray_dm)
gg_bray_pcoa <- plot_ordination(parkinsons_final_depression, pcoa_bray, color = "depression_binned_Disease") +
  labs(col = "Depression and Disease Status") + theme_bw() + stat_ellipse(level = 0.95) +
  ggtitle("Bray Curtis") + theme(plot.title = element_text(hjust = 0.5)) +
  ggside::geom_xsideboxplot(aes(fill = depression_binned_Disease, y = depression_binned_Disease), 
                            orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = depression_binned_Disease, x = depression_binned_Disease), 
                            orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void() 
gg_bray_pcoa

ggsave("bray_pcoa.png"
       , gg_bray_pcoa
       , height=4, width=5)

## unweighted unifrac ##
unifrac_dm <- distance(parkinsons_final_depression, method = "unifrac")
pcoa_unifrac <- ordinate(parkinsons_final_depression, method = "PCoA", distance = unifrac_dm)
gg_unifrac_pcoa <- plot_ordination(parkinsons_final_depression, pcoa_unifrac, color = "depression_binned_Disease") +
  labs(col = "Depression and Disease Status") + theme_bw() + stat_ellipse(level = 0.95) +
  ggtitle("Unweighted Unifrac") + theme(plot.title = element_text(hjust = 0.5))+
  ggside::geom_xsideboxplot(aes(fill = depression_binned_Disease, y = depression_binned_Disease), 
                            orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = depression_binned_Disease, x = depression_binned_Disease), 
                            orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void() 
gg_unifrac_pcoa

ggsave("unifrac_pcoa.png"
       , gg_unifrac_pcoa
       , height=4, width=5)

## weighted_unifrac ##
w_unifrac_dm <- distance(parkinsons_final_depression, method ="wunifrac")
pcoa_w_unifrac <- ordinate(parkinsons_final_depression, method="PCoA", distance=w_unifrac_dm)
gg_wunifrac_pcoa <- plot_ordination(parkinsons_final_depression, pcoa_w_unifrac, color = "depression_binned_Disease") +
  labs(col = "Depression and Disease Status") + theme_bw() + stat_ellipse(level = 0.95) +
  ggtitle("Weighted Unifrac") + theme(plot.title = element_text(hjust = 0.5)) +
  ggside::geom_xsideboxplot(aes(fill = depression_binned_Disease, y = depression_binned_Disease), 
                            orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = depression_binned_Disease, x = depression_binned_Disease), 
                            orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void() 
gg_wunifrac_pcoa

ggsave("wunifrac_pcoa.png"
       , gg_wunifrac_pcoa
       , height=4, width=5)

beta_div <- grid.arrange(gg_jac_pcoa, gg_bray_pcoa, gg_unifrac_pcoa, gg_wunifrac_pcoa)
ggsave(filename = "grid_depression_disease_diversity.png", beta_div, height = 10, width = 19)
#### Taxonomy bar plots ####

# Convert to relative abundance
parkinsons_RA <- transform_sample_counts(parkinsons_final_depression, function(x) x/sum(x))

# To remove black bars, "glom" by phylum first
parkinsons_phylum <- tax_glom(parkinsons_RA, taxrank = "Phylum", NArm=FALSE)

#Taxonomy for Disease Status
gg_taxa_disease <- plot_bar(parkinsons_phylum, fill="Phylum") + 
  facet_wrap(.~Disease, scales = "free_x")
gg_taxa_disease

#Taxonomy for Depression Status
gg_taxa_depression <- plot_bar(parkinsons_phylum, fill="Phylum") + 
  facet_wrap(.~depression_binned, scales = "free_x")
gg_taxa_depression

#Taxonomy for Disease and Depression Status 
gg_taxa_dep_PD <- plot_bar(parkinsons_phylum, fill="Phylum") + 
  facet_wrap(.~depression_binned_Disease, scales = "free_x")
gg_taxa_dep_PD

ggsave("plot_taxonomy_disease_bin.png"
       , gg_taxa_disease
       , height=13, width =35)

ggsave("plot_taxonomy_depression_bin.png"
       , gg_taxa_depression
       , height=13, width =35)

ggsave("plot_taxonomy_depPD_bin.png"
       , gg_taxa_dep_PD
       , height=13, width =35)

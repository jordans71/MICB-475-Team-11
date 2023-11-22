library(phyloseq)
library(ape)
library(tidyverse)
library(picante)
library(gridExtra)
library(ggside)

### Load in RData ###
load("parkinsons_dep_disease_final.RData")
load("parkinsons_dep_disease_rare.RData")

### Alpha Diversity ###
## Depression Status Only ##
gg_richness <- plot_richness(parkinsons_rare, x = "depression_binned") + 
  xlab("Depression_Status") + geom_boxplot()
gg_richness

ggsave(filename = "plot_richness.png", gg_richness, height = 4, width = 10)

## Depression and PD Status ##
gg_richness_dep_PD <- plot_richness(parkinsons_rare, x = "depression_binned_Disease") + 
  xlab("Depression_PD_Status") + geom_boxplot()
gg_richness_dep_PD

ggsave(filename = "plot_richness_depPD.png", gg_richness_dep_PD, height = 4, width = 10)

### Phylogenetic diversity ###
#calculate Faith's phylogenetic diversity as PD
phylo_dist <- pd(t(otu_table(parkinsons_rare)), phy_tree(parkinsons_rare), 
                 include.root = F)
sample_data(parkinsons_rare)$PD <- phylo_dist$PD
plot.pd <- ggplot(sample_data(parkinsons_rare), aes(depression_binned, PD)) + geom_boxplot() + 
  xlab("Subject #") + ylab("Phylogenetic Diversity")
plot.pd

### Beta Diversity ###
## Jaccard ##
jac_dm <- distance(parkinsons_rare, method = "jaccard", binary = T)
pcoa_jac <- ordinate(parkinsons_rare, method = "PCoa", distance = jac_dm)
gg_jac_pcoa <- plot_ordination(parkinsons_rare, pcoa_jac, color = "depression_binned") +
  labs(col = "Depression Status") + theme_bw() + stat_ellipse(level = 0.95) 
gg_jac_pcoa

ggsave("jaccard_pcoa.png"
       , gg_jac_pcoa
       , height=4, width=5)

## bray curtis ##
bray_dm <- distance(parkinsons_rare, method = "bray")
pcoa_bray <- ordinate(parkinsons_rare, method = "PCoA", distance = bray_dm)
gg_bray_pcoa <- plot_ordination(parkinsons_rare, pcoa_bray, color = "depression_binned") +
  labs(col = "Depression Status") + theme_bw() + stat_ellipse(level = 0.95) +
  ggtitle("Bray Curtis") + theme(plot.title = element_text(hjust = 0.5)) +
  ggside::geom_xsideboxplot(aes(fill = depression_binned, y = depression_binned), 
                            orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = depression_binned, x = depression_binned), 
                            orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void() 
gg_bray_pcoa

ggsave("bray_pcoa.png"
       , gg_bray_pcoa
       , height=4, width=5)

## unweighted unifrac ##
unifrac_dm <- distance(parkinsons_rare, method = "unifrac")
pcoa_unifrac <- ordinate(parkinsons_rare, method = "PCoA", distance = unifrac_dm)
gg_unifrac_pcoa <- plot_ordination(parkinsons_rare, pcoa_unifrac, color = "depression_binned") +
  labs(col = "Depression Status") + theme_bw() + stat_ellipse(level = 0.95) +
  labs(col = "Depression Status") + theme_bw() + stat_ellipse(level = 0.95) +
  ggtitle("Unweighted Unifrac") + theme(plot.title = element_text(hjust = 0.5)) +
  ggside::geom_xsideboxplot(aes(fill = depression_binned, y = depression_binned), 
                            orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = depression_binned, x = depression_binned), 
                            orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void() 
gg_unifrac_pcoa

ggsave("unifrac_pcoa.png"
       , gg_unifrac_pcoa
       , height=4, width=5)

## weighted_unifrac ##
w_unifrac_dm <- distance(parkinsons_rare, method ="wunifrac")
pcoa_w_unifrac <- ordinate(parkinsons_rare, method="PCoA", distance=w_unifrac_dm)
gg_wunifrac_pcoa <- plot_ordination(parkinsons_rare, pcoa_w_unifrac, color = "depression_binned") +
  labs(col = "Depression Status") + theme_bw() + stat_ellipse(level = 0.95) +
  labs(col = "Depression Status") + theme_bw() + stat_ellipse(level = 0.95) +
  ggtitle("Weighted Unifrac") + theme(plot.title = element_text(hjust = 0.5)) +
  ggside::geom_xsideboxplot(aes(fill = depression_binned, y = depression_binned), 
                            orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = depression_binned, x = depression_binned), 
                            orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void() 
gg_wunifrac_pcoa

ggsave("wunifrac_pcoa.png"
       , gg_wunifrac_pcoa
       , height=4, width=5)

beta_div <- grid.arrange(gg_bray_pcoa, gg_unifrac_pcoa, gg_wunifrac_pcoa)

#### Taxonomy bar plots ####

# Convert to relative abundance
parkinsons_RA <- transform_sample_counts(parkinsons_rare, function(x) x/sum(x))

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

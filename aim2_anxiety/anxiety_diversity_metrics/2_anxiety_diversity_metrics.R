library(phyloseq)
library(ape)
library(tidyverse)
library(picante)
library(vegan)
library(gridExtra)

#### Load in RData ####
load("parkinsons_final_anxiety.RData")
load("parkinsons2_filt_nolow_samps.RData")
load("parkinsons2_rare.RData")
load("parkinsons_final_depression.RData")
load("parkinsons2.RData")

### Alpha Diversity of PD Patients ###
#filter for PD patients and remove NAs
PD_patients <- subset_samples(parkinsons_rare_new_anxiety, `Disease` == "PD", !is.na(anxitey_binned))

#Plot Alpha Diversity Metrics for PD patients with anxiety 
gg_richness <- plot_richness(PD_patients, x = "anxitey_binned") +
  xlab("PD_anxiety_level") +geom_boxplot() + ggtitle("PD Anxiety Alpha Diversity Metrics ") + 
  theme(plot.title = element_text(hjust = 0.5)) + xlab("PD Patient Anxiety Level") 
gg_richness

ggsave(filename = "Alpha_diversity_PD_anxiety_level_new.png"
       , gg_richness
       , height=4, width=6)

#Alpha_Diversity with control patients 
#filter for control patients and remove NAs
ctrl_patients <- subset_samples(parkinsons_rare_new_anxiety, `Disease` == "Control", !is.na(anxitey_binned))

gg_richness_control <- plot_richness(ctrl_patients, x = "anxitey_binned") +
  xlab("Control_anxiety_level") + geom_boxplot() + ggtitle("Control Anxiety Alpha Diversity Metrics") +
  theme(plot.title = element_text(hjust = 0.5)) + xlab("Healthy Control Anxiety Level")
gg_richness_control

alpha_anxiety <- grid.arrange(gg_richness, gg_richness_control, ncol =1)
ggsave(filename = "alpha_anxiety.png", alpha_anxiety)
### Beta Diversity ###
## Jaccard ## 
#PD patients
jac_dm <- distance(PD_patients, method = "jaccard", binary = T)
pcoa_jac_PD <- ordinate(PD_patients, method = "NMDS", distance = jac_dm)
PD_anxiety_jac <- plot_ordination(PD_patients, pcoa_jac_PD, color = "anxitey_binned") +
  labs(col = "Anxiety Level") + theme_bw() + stat_ellipse(level = 0.95) +
  ggside::geom_xsideboxplot(aes(fill = anxitey_binned, y = anxitey_binned), orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = anxitey_binned, x = anxitey_binned), orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void() 
PD_anxiety_jac

#Healthy Controls
jac_dm_ctrl <- distance(ctrl_patients, method = "jaccard", binary = T)
pcoa_jac_ctrl <- ordinate(ctrl_patients, method = "NMDS", distance = jac_dm_ctrl)
ctrl_anxiety_jac <- plot_ordination(ctrl_patients, pcoa_jac_ctrl, color = "anxitey_binned") + 
  labs(col = "Anxiety Level") + theme_bw() + stat_ellipse(level = 0.95) +
  ggside::geom_xsideboxplot(aes(fill = anxitey_binned, y = anxitey_binned), orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = anxitey_binned, x = anxitey_binned), orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void() 
ctrl_anxiety_jac

##Bray-Curtis ##
# PD patients
bc_dm <- distance(PD_patients, method="bray")
pcoa_bc <- ordinate(PD_patients, method="PCoA", distance=bc_dm)
PD_anxiety_bray <- plot_ordination(PD_patients, pcoa_bc, color = "anxitey_binned") + 
  labs(col = "anxiety level") + theme_bw() + stat_ellipse(level = 0.95) +
  ggside::geom_xsideboxplot(aes(fill = anxitey_binned, y = anxitey_binned), orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = anxitey_binned, x = anxitey_binned), orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void()
PD_anxiety_bray

## Healthy controls
bc_dm_ctrl <- distance(ctrl_patients, method="bray")
pcoa_bc_ctrl <- ordinate(ctrl_patients, method="PCoA", distance=bc_dm)
ctrl_anxiety_bray <- plot_ordination(ctrl_patients, pcoa_bc_ctrl, color = "anxitey_binned") + 
  labs(col = "anxiety level") + theme_bw() + stat_ellipse(level = 0.95) +
  ggside::geom_xsideboxplot(aes(fill = anxitey_binned, y = anxitey_binned), orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = anxitey_binned, x = anxitey_binned), orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void()
ctrl_anxiety_beta_diversity

## Unweighted Unifrac ##
#PD patients 
unifrac_dm_PD <- distance(PD_patients, method = "unifrac")
pcoa_unifrac_PD <- ordinate(PD_patients, method = "PCoA", distance = unifrac_dm_PD)
gg_unifrac_pcoa <- plot_ordination(PD_patients, pcoa_unifrac_PD, color = "anxitey_binned") +
  labs(col = "Anxiety Level") + theme_bw() + stat_ellipse(level = 0.95) +
  ggtitle("Unweighted Unifrac") + theme(plot.title = element_text(hjust = 0.5)) +
  ggside::geom_xsideboxplot(aes(fill = anxitey_binned, y = anxitey_binned), orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = anxitey_binned, x = anxitey_binned), orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void() 
gg_unifrac_pcoa

#Healthy Controls
unifrac_dm_ctrl <- distance(ctrl_patients, method = "unifrac")
pcoa_unifrac_ctrl <- ordinate(ctrl_patients, method = "PCoA", distance = unifrac_dm_ctrl)
gg_unifrac_pcoa_ctrl <- plot_ordination(ctrl_patients, pcoa_unifrac_ctrl, color = "anxitey_binned") +
  labs(col = "Anxiety Level") + theme_bw() + stat_ellipse(level = 0.95) +
  ggtitle("Unweighted Unifrac") + theme(plot.title = element_text(hjust = 0.5)) +
  ggside::geom_xsideboxplot(aes(fill = anxitey_binned, y = anxitey_binned), orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = anxitey_binned, x = anxitey_binned), orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void() 
gg_unifrac_pcoa_ctrl

## Weighted Unifrac ##
#PD patients 
w_unifrac_dm_PD <- distance(PD_patients, method ="wunifrac")
pcoa_w_unifrac_PD <- ordinate(PD_patients, method="PCoA", distance=w_unifrac_dm_PD)
gg_wunifrac_pcoa_PD <- plot_ordination(PD_patients, pcoa_w_unifrac_PD, color = "anxitey_binned") +
  labs(col = "Anxiety Status") + theme_bw() + stat_ellipse(level = 0.95) +
  ggtitle("Weighted Unifrac") + theme(plot.title = element_text(hjust = 0.5)) +
  ggside::geom_xsideboxplot(aes(fill = anxitey_binned, y = anxitey_binned), orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = anxitey_binned, x = anxitey_binned), orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void() 
gg_wunifrac_pcoa_PD

#Healthy Controls
w_unifrac_dm_ctrl <- distance(ctrl_patients, method ="wunifrac")
pcoa_w_unifrac_ctrl <- ordinate(ctrl_patients, method="PCoA", distance=w_unifrac_dm_ctrl)
gg_wunifrac_pcoa_ctrl <- plot_ordination(ctrl_patients, pcoa_w_unifrac_ctrl, color = "anxitey_binned") +
  labs(col = "Anxiety Status") + theme_bw() + stat_ellipse(level = 0.95) +
  ggtitle("Weighted Unifrac") + theme(plot.title = element_text(hjust = 0.5)) +
  ggside::geom_xsideboxplot(aes(fill = anxitey_binned, y = anxitey_binned), orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = anxitey_binned, x = anxitey_binned), orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void() 
gg_wunifrac_pcoa_ctrl

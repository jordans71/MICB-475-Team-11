library(phyloseq)
library(ape)
library(tidyverse)
library(picante)
library(vegan)
library(gridExtra)

#### Load in RData ####
load("parkinsons_final_sleep.RData")

### Alpha Diversity of PD Patients ###
#filter for PD patients and remove NAs
PD_patients <- subset_samples(parkinsons_final_sleep, `Disease` == "PD", !is.na(Sleep_problems))

#Plot Alpha Diversity Metrics for PD patients with sleep problems
gg_richness_sleep <- plot_richness(PD_patients, x = "Sleep_problems") +
  xlab("PD Sleep Problems") +geom_boxplot() + ggtitle("PD Sleep Alpha Diversity Metrics ") + 
  theme(plot.title = element_text(hjust = 0.5))  
gg_richness_sleep

ggsave(filename = "Alpha_diversity_PD_sleep.png"
       , gg_richness_sleep)

#Alpha_Diversity with control patients 
#filter for control patients and remove NAs
ctrl_patients <- subset_samples(parkinsons_final_sleep, `Disease` == "Control", !is.na(Sleep_problems))

gg_richness_ctrl_sleep <- plot_richness(ctrl_patients, x = "Sleep_problems") +
  xlab("Control Sleep Problems") + geom_boxplot() + ggtitle("Control Sleep Alpha Diversity Metrics") +
  theme(plot.title = element_text(hjust = 0.5)) 
gg_richness_ctrl_sleep

ggsave(filename = "Alpha_diversity_ctrl_sleep.png"
       , gg_richness_ctrl_sleep)

alpha_sleep <- grid.arrange(gg_richness_sleep, gg_richness_ctrl_sleep, ncol = 1)
ggsave(filename = "alpha_sleep_grid.png", alpha_sleep)

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

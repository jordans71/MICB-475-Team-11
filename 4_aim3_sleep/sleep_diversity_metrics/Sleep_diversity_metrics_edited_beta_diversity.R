### Load libraries ###
library(phyloseq)
library(ape) 
library(tidyverse)
library(vegan)
### Load sleep_problems  phyloseq Object 
load("parkinsons_final_sleep_problems.RData")



#Filter out PD patients 
PD_patients <- subset_samples(parkinsons_final_sleep_problems , `Disease` == "PD")

#Filter out control
Ctrl_patients <- subset_samples(parkinsons_final_sleep_problems , `Disease` == "Control")

###Data frame 
#PD
samp_dat_wdiv_PD <- data.frame(sample_data(PD_patients), estimate_richness(PD_patients))

#Control
samp_dat_wdiv_Ctrl <- data.frame(sample_data(Ctrl_patients), estimate_richness(Ctrl_patients))



### Beta Diversity ###
## Jaccard ## 
#PD patients
jac_dm <- distance(PD_patients, method = "jaccard", binary = T)
pcoa_jac_PD <- ordinate(PD_patients, method = "NMDS", distance = jac_dm)
PD_sleep_problems_jac <- plot_ordination(PD_patients, pcoa_jac_PD, color = "Sleep_problems") +
  labs(col = "sleep_problems_Level") + theme_bw() + stat_ellipse(level = 0.95) +
  ggtitle("Jaccard_PD_sleep_problems") + theme(plot.title = element_text(hjust = 0.5)) +
  ggside::geom_xsideboxplot(aes(fill = Sleep_problems , y = Sleep_problems), orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = Sleep_problems , x = Sleep_problems), orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void() 
PD_sleep_problems_jac
ggsave("PD_sleep_problems_jac_pcoa.png"
       , PD_sleep_problems_jac
       , height=4, width=6)

adonis2(jac_dm ~ `Sleep_problems`, data = samp_dat_wdiv_PD)


#Healthy Controls
jac_dm_ctrl <- distance(Ctrl_patients, method = "jaccard", binary = T)
pcoa_jac_ctrl <- ordinate(Ctrl_patients, method = "NMDS", distance = jac_dm_ctrl)
ctrl_sleep_problems_jac <- plot_ordination(Ctrl_patients, pcoa_jac_ctrl, color = "Sleep_problems") + 
  labs(col = "sleep_problems_Level") + theme_bw() + stat_ellipse(level = 0.95) +
  ggtitle("Jaccard_Control_sleep_problems") + theme(plot.title = element_text(hjust = 0.5)) +
  ggside::geom_xsideboxplot(aes(fill = Sleep_problems , y = Sleep_problems), orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = Sleep_problems , x = Sleep_problems), orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void() 
ctrl_sleep_problems_jac

ggsave("ctrl_sleep_problems _jac_pcoa.png"
       , ctrl_sleep_problems_jac
       , height=4, width=6)
adonis2(jac_dm_ctrl ~ `Sleep_problems`, data = samp_dat_wdiv_Ctrl)

##Bray-Curtis ##
# PD patients
bc_dm <- distance(PD_patients, method="bray")
pcoa_bc_PD <- ordinate(PD_patients, method="PCoA", distance=bc_dm)
PD_sleep_problems_bray <- plot_ordination(PD_patients, pcoa_bc_PD, color = "Sleep_problems") + 
  labs(col = "sleep_problems_level") + theme_bw() + stat_ellipse(level = 0.95) +
  ggtitle("Bray_Curtis_PD_sleep_problems") + theme(plot.title = element_text(hjust = 0.5)) +
  ggside::geom_xsideboxplot(aes(fill = Sleep_problems , y = Sleep_problems ), orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = Sleep_problems , x = Sleep_problems ), orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void()
PD_sleep_problems_bray
ggsave("PD_sleep_problems_bray_pcoa.png"
       , PD_sleep_problems_bray
       , height=4, width=6)

adonis2(bc_dm ~ `Sleep_problems`, data = samp_dat_wdiv_PD)


## Healthy controls
bc_dm_ctrl <- distance(Ctrl_patients, method="bray")
pcoa_bc_ctrl <- ordinate(Ctrl_patients, method="PCoA", distance=bc_dm_ctrl)
ctrl_sleep_problems_bray <- plot_ordination(Ctrl_patients, pcoa_bc_ctrl, color = "Sleep_problems") + 
  labs(col = "sleep_problems_level") + theme_bw() + stat_ellipse(level = 0.95) +
  ggtitle("Bray_Curtis_Control_sleep_problems") + theme(plot.title = element_text(hjust = 0.5)) +
  ggside::geom_xsideboxplot(aes(fill = Sleep_problems , y = Sleep_problems ), orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = Sleep_problems , x = Sleep_problems ), orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void()
ctrl_sleep_problems_bray
ggsave("ctrl_sleep_problems_bray_pcoa.png"
       , ctrl_sleep_problems_bray
       , height=4, width=6)
adonis2(bc_dm_ctrl ~ `Sleep_problems`, data = samp_dat_wdiv_Ctrl)

## Unweighted Unifrac ##
#PD patients 
unifrac_dm_PD <- distance(PD_patients, method = "unifrac")
pcoa_unifrac_PD <- ordinate(PD_patients, method = "PCoA", distance = unifrac_dm_PD)
PD_gg_unifrac_pcoa <- plot_ordination(PD_patients, pcoa_unifrac_PD, color = "Sleep_problems") +
  labs(col = "sleep_problems_Level") + theme_bw() + stat_ellipse(level = 0.95) +
  ggtitle("Unweighte_Unifrac_PD_sleep_problems") + theme(plot.title = element_text(hjust = 0.5)) +
  ggside::geom_xsideboxplot(aes(fill = Sleep_problems , y = Sleep_problems ), orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = Sleep_problems , x = Sleep_problems ), orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void() 
PD_gg_unifrac_pcoa
ggsave("PD_sleep_problems _unifrac.png"
       , PD_gg_unifrac_pcoa
       , height=4, width=6)

adonis2(unifrac_dm_PD ~ `Sleep_problems`, data = samp_dat_wdiv_PD)

#Healthy Controls
unifrac_dm_ctrl <- distance(Ctrl_patients, method = "unifrac")
pcoa_unifrac_ctrl <- ordinate(Ctrl_patients, method = "PCoA", distance = unifrac_dm_ctrl)
ctrl_gg_unifrac_pcoa <- plot_ordination(Ctrl_patients, pcoa_unifrac_ctrl, color = "Sleep_problems") +
  labs(col = "sleep_problems_Level") + theme_bw() + stat_ellipse(level = 0.95) +
  ggtitle("Unweighted_Unifrac_Control_sleep_problems") + theme(plot.title = element_text(hjust = 0.5)) +
  ggside::geom_xsideboxplot(aes(fill = Sleep_problems , y = Sleep_problems), orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = Sleep_problems , x = Sleep_problems), orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void() 
ctrl_gg_unifrac_pcoa
ggsave("ctrl_sleep_problems_unifrac_pcoa.png"
       , ctrl_gg_unifrac_pcoa
       , height=4, width=6)

adonis2(unifrac_dm_ctrl ~ `Sleep_problems`, data = samp_dat_wdiv_Ctrl)

## Weighted Unifrac ##
#PD patients 
w_unifrac_dm_PD <- distance(PD_patients, method ="wunifrac")
pcoa_w_unifrac_PD <- ordinate(PD_patients, method="PCoA", distance=w_unifrac_dm_PD)
gg_wunifrac_pcoa_PD <- plot_ordination(PD_patients, pcoa_w_unifrac_PD, color = "Sleep_problems") +
  labs(col = "sleep_problems_level") + theme_bw() + stat_ellipse(level = 0.95) +
  ggtitle("Weighted_Unifrac_PD_sleep_problems") + theme(plot.title = element_text(hjust = 0.5)) +
  ggside::geom_xsideboxplot(aes(fill = Sleep_problems , y = Sleep_problems ), orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = Sleep_problems , x = Sleep_problems ), orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void() 
gg_wunifrac_pcoa_PD

ggsave("PD_sleep_problems _wunifrac_pcoa.png"
       , gg_wunifrac_pcoa_PD
       , height=4, width=6)

adonis2(w_unifrac_dm_PD ~ `Sleep_problems`, data = samp_dat_wdiv_PD)


#Healthy Controls
w_unifrac_dm_ctrl <- distance(Ctrl_patients, method ="wunifrac")
pcoa_w_unifrac_ctrl <- ordinate(Ctrl_patients, method="PCoA", distance=w_unifrac_dm_ctrl)
gg_wunifrac_pcoa_ctrl <- plot_ordination(Ctrl_patients, pcoa_w_unifrac_ctrl, color = "Sleep_problems") +
  labs(col = "sleep_problems_level") + theme_bw() + stat_ellipse(level = 0.95) +
  ggtitle("Weighted_Unifrac_Control_sleep_problems") + theme(plot.title = element_text(hjust = 0.5)) +
  ggside::geom_xsideboxplot(aes(fill = Sleep_problems , y = Sleep_problems), orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = Sleep_problems , x = Sleep_problems), orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void() 
gg_wunifrac_pcoa_ctrl
ggsave("ctrl_sleep_problems_wunifrac_pcoa.png"
       , gg_wunifrac_pcoa_ctrl
       , height=4, width=6)

adonis2(w_unifrac_dm_ctrl ~ `Sleep_problems`, data = samp_dat_wdiv_Ctrl)


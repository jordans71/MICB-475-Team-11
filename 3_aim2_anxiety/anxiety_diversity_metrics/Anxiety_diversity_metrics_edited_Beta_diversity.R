### Load libraries ###
library(phyloseq)
library(ape) 
library(tidyverse)
library(vegan)
### Load anxiety phyloseq Object 
load("parkinsons_final_anxiety.RData")



#Filter out PD patients 
PD_patients <- subset_samples(parkinsons_final_anxiety, `Disease` == "PD")

#Filter out control
Ctrl_patients <- subset_samples(parkinsons_final_anxiety, `Disease` == "Control")

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
PD_anxiety_jac <- plot_ordination(PD_patients, pcoa_jac_PD, color = "anxiety_binned") +
  labs(col = "Anxiety Level") + theme_bw() + stat_ellipse(level = 0.95) +
  ggtitle("Jaccard_PD_anxiety") + theme(plot.title = element_text(hjust = 0.5)) +
  ggside::geom_xsideboxplot(aes(fill = anxiety_binned, y = anxiety_binned), orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = anxiety_binned, x = anxiety_binned), orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void() 
PD_anxiety_jac
ggsave("PD_anxiety_jac_pcoa.png"
       , PD_anxiety_jac
       , height=4, width=6)

adonis2(jac_dm ~ `anxiety_binned`, data = samp_dat_wdiv_PD)


#Healthy Controls
jac_dm_ctrl <- distance(Ctrl_patients, method = "jaccard", binary = T)
pcoa_jac_ctrl <- ordinate(Ctrl_patients, method = "NMDS", distance = jac_dm_ctrl)
ctrl_anxiety_jac <- plot_ordination(Ctrl_patients, pcoa_jac_ctrl, color = "anxiety_binned") + 
  labs(col = "Anxiety Level") + theme_bw() + stat_ellipse(level = 0.95) +
  ggtitle("Jaccard_Control_anxiety") + theme(plot.title = element_text(hjust = 0.5)) +
  ggside::geom_xsideboxplot(aes(fill = anxiety_binned, y = anxiety_binned), orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = anxiety_binned, x = anxiety_binned), orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void() 
ctrl_anxiety_jac

ggsave("ctrl_anxiety_jac_pcoa.png"
       , ctrl_anxiety_jac
       , height=4, width=6)
adonis2(jac_dm_ctrl ~ `anxiety_binned`, data = samp_dat_wdiv_Ctrl)

##Bray-Curtis ##
# PD patients
bc_dm <- distance(PD_patients, method="bray")
pcoa_bc_PD <- ordinate(PD_patients, method="PCoA", distance=bc_dm)
PD_anxiety_bray <- plot_ordination(PD_patients, pcoa_bc_PD, color = "anxiety_binned") + 
  labs(col = "anxiety level") + theme_bw() + stat_ellipse(level = 0.95) +
  ggtitle("Bray_Curtis_PD_anxiety") + theme(plot.title = element_text(hjust = 0.5)) +
  ggside::geom_xsideboxplot(aes(fill = anxiety_binned, y = anxiety_binned), orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = anxiety_binned, x = anxiety_binned), orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void()
PD_anxiety_bray
ggsave("PD_anxiety_bray_pcoa.png"
       , PD_anxiety_bray
       , height=4, width=6)

adonis2(bc_dm ~ `anxiety_binned`, data = samp_dat_wdiv_PD)


## Healthy controls
bc_dm_ctrl <- distance(Ctrl_patients, method="bray")
pcoa_bc_ctrl <- ordinate(Ctrl_patients, method="PCoA", distance=bc_dm_ctrl)
ctrl_anxiety_bray <- plot_ordination(Ctrl_patients, pcoa_bc_ctrl, color = "anxiety_binned") + 
  labs(col = "anxiety level") + theme_bw() + stat_ellipse(level = 0.95) +
  ggtitle("Bray_Curtis_Control_anxiety") + theme(plot.title = element_text(hjust = 0.5)) +
  ggside::geom_xsideboxplot(aes(fill = anxiety_binned, y = anxiety_binned), orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = anxiety_binned, x = anxiety_binned), orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void()
ctrl_anxiety_bray
ggsave("ctrl_anxiety_bray_pcoa.png"
       , ctrl_anxiety_bray
       , height=4, width=6)
adonis2(bc_dm_ctrl ~ `anxiety_binned`, data = samp_dat_wdiv_Ctrl)

## Unweighted Unifrac ##
#PD patients 
unifrac_dm_PD <- distance(PD_patients, method = "unifrac")
pcoa_unifrac_PD <- ordinate(PD_patients, method = "PCoA", distance = unifrac_dm_PD)
PD_gg_unifrac_pcoa <- plot_ordination(PD_patients, pcoa_unifrac_PD, color = "anxiety_binned") +
  labs(col = "Anxiety Level") + theme_bw() + stat_ellipse(level = 0.95) +
  ggtitle("Unweighted_Unifrac_PD_anxiety") + theme(plot.title = element_text(hjust = 0.5)) +
  ggside::geom_xsideboxplot(aes(fill = anxiety_binned, y = anxiety_binned), orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = anxiety_binned, x = anxiety_binned), orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void() 
PD_gg_unifrac_pcoa
ggsave("PD_anxiety_unifrac.png"
       , PD_gg_unifrac_pcoa
       , height=4, width=6)

adonis2(unifrac_dm_PD ~ `anxiety_binned`, data = samp_dat_wdiv_PD)

#Healthy Controls
unifrac_dm_ctrl <- distance(Ctrl_patients, method = "unifrac")
pcoa_unifrac_ctrl <- ordinate(Ctrl_patients, method = "PCoA", distance = unifrac_dm_ctrl)
ctrl_gg_unifrac_pcoa <- plot_ordination(Ctrl_patients, pcoa_unifrac_ctrl, color = "anxiety_binned") +
  labs(col = "Anxiety Level") + theme_bw() + stat_ellipse(level = 0.95) +
  ggtitle("Unweighted_Unifrac_Control_anxiety") + theme(plot.title = element_text(hjust = 0.5)) +
  ggside::geom_xsideboxplot(aes(fill = anxiety_binned, y = anxiety_binned), orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = anxiety_binned, x = anxiety_binned), orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void() 
ctrl_gg_unifrac_pcoa
ggsave("ctrl_anxiety_unifrac_pcoa.png"
       , ctrl_gg_unifrac_pcoa
       , height=4, width=6)

adonis2(unifrac_dm_ctrl ~ `anxiety_binned`, data = samp_dat_wdiv_Ctrl)

## Weighted Unifrac ##
#PD patients 
w_unifrac_dm_PD <- distance(PD_patients, method ="wunifrac")
pcoa_w_unifrac_PD <- ordinate(PD_patients, method="PCoA", distance=w_unifrac_dm_PD)
gg_wunifrac_pcoa_PD <- plot_ordination(PD_patients, pcoa_w_unifrac_PD, color = "anxiety_binned") +
  labs(col = "Anxiety Status") + theme_bw() + stat_ellipse(level = 0.95) +
  ggtitle("Weighted_Unifrac_PD_anxiety") + theme(plot.title = element_text(hjust = 0.5)) +
  ggside::geom_xsideboxplot(aes(fill = anxiety_binned, y = anxiety_binned), orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = anxiety_binned, x = anxiety_binned), orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void() 
gg_wunifrac_pcoa_PD

ggsave("PD_anxiety_wunifrac_pcoa.png"
       , gg_wunifrac_pcoa_PD
       , height=4, width=6)

adonis2(w_unifrac_dm_PD ~ `anxiety_binned`, data = samp_dat_wdiv_PD)


#Healthy Controls
w_unifrac_dm_ctrl <- distance(Ctrl_patients, method ="wunifrac")
pcoa_w_unifrac_ctrl <- ordinate(Ctrl_patients, method="PCoA", distance=w_unifrac_dm_ctrl)
gg_wunifrac_pcoa_ctrl <- plot_ordination(Ctrl_patients, pcoa_w_unifrac_ctrl, color = "anxiety_binned") +
  labs(col = "Anxiety Status") + theme_bw() + stat_ellipse(level = 0.95) +
  ggtitle("Weighted_Unifrac_Control_anxiety") + theme(plot.title = element_text(hjust = 0.5)) +
  ggside::geom_xsideboxplot(aes(fill = anxiety_binned, y = anxiety_binned), orientation = "y") +
  ggside::geom_ysideboxplot(aes(fill = anxiety_binned, x = anxiety_binned), orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void() 
gg_wunifrac_pcoa_ctrl
ggsave("ctrl_anxiety_wunifrac_pcoa.png"
       , gg_wunifrac_pcoa_ctrl
       , height=4, width=6)

adonis2(w_unifrac_dm_ctrl ~ `anxiety_binned`, data = samp_dat_wdiv_Ctrl)


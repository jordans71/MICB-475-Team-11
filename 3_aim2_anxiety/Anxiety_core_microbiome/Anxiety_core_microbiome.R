library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)
library(sf)

#load data
load("parkinsons2_filt_nolow_samps.RData")

# make phyloseq object for PD and one for Control, then make 6 total venn diagrams
# 1 for each aim using PD and 1 for Control

# Keep only Control
parkinsons_control <- subset_samples(parkinsons_filt_nolow_samps, `Disease`=="Control")
View(sample_data(parkinsons_control))
save(parkinsons_control, file="parkinsons_control.RData")

# Keep only PD
parkinsons_PD <- subset_samples(parkinsons_filt_nolow_samps, `Disease` == "PD")
View(sample_data(parkinsons_PD))
save(parkinsons_PD, file="parkinsons_PD.RData")


# Convert to relative abundance
parkinsons_RA <- transform_sample_counts(parkinsons_filt_nolow_samps, fun=function(x) x/sum(x))
parkinsons_control_RA <- transform_sample_counts(parkinsons_control, fun=function(x) x/sum(x))
parkinsons_PD_RA <- transform_sample_counts(parkinsons_PD, fun=function(x) x/sum(x))
View(sample_data(parkinsons_PD_RA))


##### ANXIETY CONTROL #####

# Filter dataset by anxiety
parkinsons_control_hianx <- subset_samples(parkinsons_control_RA, `anxitey_binned`=="High")
parkinsons_control_loanx <- subset_samples(parkinsons_control_RA, `anxitey_binned`=="Low")

# What ASVs are found in more than 25% of samples in each anxiety category?

control_hianx_ASVs <- core_members(parkinsons_control_hianx, detection=0.01, prevalence = 0.25)
control_loanx_ASVs <- core_members(parkinsons_control_loanx, detection=0.01, prevalence = 0.25)

# What are these ASVs?
prune_taxa(control_hianx_ASVs,parkinsons_filt_nolow_samps) %>%
  tax_table()

# can plot those ASVs' relative abundance
prune_taxa(control_loanx_ASVs,parkinsons_control_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`anxitey_binned`, scales ="free")

# What if we wanted to make a Venn diagram of all the ASVs that showed up in each treatment?
control_hianx_list <- core_members(parkinsons_control_hianx, detection=0.01, prevalence = 0.25)
control_loanx_list <- core_members(parkinsons_control_loanx, detection=0.01, prevalence = 0.25)

control_anx_list_full <- list(High_Anxiety = control_hianx_list, Low_Anxiety = control_loanx_list)

# Create a Venn diagram using all the ASVs shared and unique to anxious and non-anxious groups
control_anxiety_venn <- ggVennDiagram(x = control_anx_list_full)
control_anxiety_venn

ggsave("control_anxiety_venn.png", control_anxiety_venn)


##### ANXIETY PD #####

# Filter dataset by anxiety
parkinsons_PD_hianx <- subset_samples(parkinsons_PD_RA, `anxitey_binned`=="High")
parkinsons_PD_loanx <- subset_samples(parkinsons_PD_RA, `anxitey_binned`=="Low")

# What ASVs are found in more than 25% of samples in each anxiety category?

PD_hianx_ASVs <- core_members(parkinsons_PD_hianx, detection=0.01, prevalence = 0.25)
PD_loanx_ASVs <- core_members(parkinsons_PD_loanx, detection=0.01, prevalence = 0.25)

# What are these ASVs?
prune_taxa(PD_hianx_ASVs,parkinsons_filt_nolow_samps) %>%
  tax_table()

# can plot those ASVs' relative abundance
prune_taxa(PD_loanx_ASVs,parkinsons_PD_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`anxitey_binned`, scales ="free")

# What if we wanted to make a Venn diagram of all the ASVs that showed up in each treatment?
PD_hianx_list <- core_members(parkinsons_PD_hianx, detection=0.01, prevalence = 0.25)
PD_loanx_list <- core_members(parkinsons_PD_loanx, detection=0.01, prevalence = 0.25)

PD_anx_list_full <- list(High_Anxiety = PD_hianx_list, Low_Anxiety = PD_loanx_list)

# Create a Venn diagram using all the ASVs shared and unique to depressed and non-depressed groups
PD_anxiety_venn <- ggVennDiagram(x = PD_anx_list_full)
PD_anxiety_venn

ggsave("PD_anxiety_venn.png", PD_anxiety_venn)


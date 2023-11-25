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

##### DEPRESSION CONTROL #####

# Filter dataset by depression
parkinsons_control_hidep <- subset_samples(parkinsons_control_RA, `depression_binned`=="High")
parkinsons_control_lodep <- subset_samples(parkinsons_control_RA, `depression_binned`=="Low")

# What ASVs are found in more than 25% of samples in each depression category?

control_hidep_ASVs <- core_members(parkinsons_control_hidep, detection=0.01, prevalence = 0.25)
control_lodep_ASVs <- core_members(parkinsons_control_lodep, detection=0.01, prevalence = 0.25)

# What are these ASVs?
prune_taxa(control_hidep_ASVs,parkinsons_filt_nolow_samps) %>%
  tax_table()

# can plot those ASVs' relative abundance
prune_taxa(control_lodep_ASVs,parkinsons_control_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`depression_binned`, scales ="free")

# What if we wanted to make a Venn diagram of all the ASVs that showed up in each treatment?
control_hidep_list <- core_members(parkinsons_control_hidep, detection=0.01, prevalence = 0.25)
control_lodep_list <- core_members(parkinsons_control_lodep, detection=0.01, prevalence = 0.25)

control_dep_list_full <- list(High_Depresssion = control_hidep_list, Low_Depression = control_lodep_list)

# Create a Venn diagram using all the ASVs shared and unique to depressed and non-depressed groups
control_depression_venn <- ggVennDiagram(x = control_dep_list_full)
control_depression_venn

ggsave("control_depression_venn.png", control_depression_venn)


##### DEPRESSION PD #####

# Filter dataset by depression
parkinsons_PD_hidep <- subset_samples(parkinsons_PD_RA, `depression_binned`=="High")
parkinsons_PD_lodep <- subset_samples(parkinsons_PD_RA, `depression_binned`=="Low")

# What ASVs are found in more than 25% of samples in each depression category?

PD_hidep_ASVs <- core_members(parkinsons_PD_hidep, detection=0.01, prevalence = 0.25)
PD_lodep_ASVs <- core_members(parkinsons_PD_lodep, detection=0.01, prevalence = 0.25)

# What are these ASVs?
prune_taxa(PD_hidep_ASVs,parkinsons_filt_nolow_samps) %>%
  tax_table()

# can plot those ASVs' relative abundance
prune_taxa(PD_lodep_ASVs,parkinsons_PD_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`depression_binned`, scales ="free")

# What if we wanted to make a Venn diagram of all the ASVs that showed up in each treatment?
PD_hidep_list <- core_members(parkinsons_PD_hidep, detection=0.01, prevalence = 0.25)
PD_lodep_list <- core_members(parkinsons_PD_lodep, detection=0.01, prevalence = 0.25)

PD_dep_list_full <- list(High_Depresssion = PD_hidep_list, Low_Depression = PD_lodep_list)

# Create a Venn diagram using all the ASVs shared and unique to depressed and non-depressed groups
PD_depression_venn <- ggVennDiagram(x = PD_dep_list_full)
PD_depression_venn

ggsave("PD_depression_venn.png", PD_depression_venn)


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



##### SLEEP PROBLEMS CONTROL #####

# Filter dataset by sleep problems
parkinsons_control_sp <- subset_samples(parkinsons_control_RA, `Sleep_problems`=="Yes")
parkinsons_control_nosp <- subset_samples(parkinsons_control_RA, `Sleep_problems`=="No")

# What ASVs are found in more than 25% of samples in each sleep category?

control_sp_ASVs <- core_members(parkinsons_control_sp, detection=0.01, prevalence = 0.25)
control_nosp_ASVs <- core_members(parkinsons_control_nosp, detection=0.01, prevalence = 0.25)

# What are these ASVs?
prune_taxa(control_sp_ASVs,parkinsons_filt_nolow_samps) %>%
  tax_table()

# can plot those ASVs' relative abundance
prune_taxa(control_nosp_ASVs,parkinsons_control_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`Sleep_problems`, scales ="free")

# What if we wanted to make a Venn diagram of all the ASVs that showed up in each treatment?
control_sp_list <- core_members(parkinsons_control_sp, detection=0.01, prevalence = 0.25)
control_nosp_list <- core_members(parkinsons_control_nosp, detection=0.01, prevalence = 0.25)

control_sp_list_full <- list(Sleep_Problems = control_sp_list, No_Sleep_Problems = control_nosp_list)

# Create a Venn diagram using all the ASVs shared and unique to sleep problem and non-sleep problem groups
control_sleep_problem_venn <- ggVennDiagram(x = control_sp_list_full)
control_sleep_problem_venn

ggsave("control_sleep_problem_venn.png", control_sleep_problem_venn)




##### SLEEP PROBLEMS PD #####

# Filter dataset by sleep problems
parkinsons_PD_sp <- subset_samples(parkinsons_PD_RA, `Sleep_problems`=="Yes")
parkinsons_PD_nosp <- subset_samples(parkinsons_PD_RA, `Sleep_problems`=="No")

# What ASVs are found in more than 25% of samples in each sleep category?

PD_sp_ASVs <- core_members(parkinsons_PD_sp, detection=0.01, prevalence = 0.25)
PD_nosp_ASVs <- core_members(parkinsons_PD_nosp, detection=0.01, prevalence = 0.25)

# What are these ASVs?
prune_taxa(PD_sp_ASVs,parkinsons_filt_nolow_samps) %>%
  tax_table()

# can plot those ASVs' relative abundance
prune_taxa(PD_nosp_ASVs,parkinsons_PD_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`Sleep_problems`, scales ="free")

# What if we wanted to make a Venn diagram of all the ASVs that showed up in each treatment?
PD_sp_list <- core_members(parkinsons_PD_sp, detection=0.01, prevalence = 0.25)
PD_nosp_list <- core_members(parkinsons_PD_nosp, detection=0.01, prevalence = 0.25)

PD_sp_list_full <- list(Sleep_Problems = PD_sp_list, No_Sleep_Problems = PD_nosp_list)

# Create a Venn diagram using all the ASVs shared and unique to sleep problem and non-sleep problem groups
PD_sleep_problem_venn <- ggVennDiagram(x = PD_sp_list_full)
PD_sleep_problem_venn

ggsave("PD_sleep_problem_venn.png", PD_sleep_problem_venn)

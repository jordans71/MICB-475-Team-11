#!/usr/bin/env Rscript
install.packages('cowplot')
install.packages('gridExtra')

library(tidyverse)
library(phyloseq)
library(DESeq2)
library(cowplot)
library(ggplot2)
library(gridExtra)


#### Load data ####
load("~/Desktop/Project_2_Parkinsons/Trial_Folder_Rae/1_make_phyloseq_objects/parkinsons_final_depression.RData")
load("~/Desktop/Project_2_Parkinsons/Trial_Folder_Rae/1_make_phyloseq_objects/parkinsons_final_anxiety.RData")
load("~/Desktop/Project_2_Parkinsons/Trial_Folder_Rae/1_make_phyloseq_objects/parkinsons_final_sleep.RData")

#### DESeq Volcano Plot for Depression ####

depression_plus1 <- transform_sample_counts(parkinsons_final_depression, function(x) x+1)
depression_deseq <- phyloseq_to_deseq2(depression_plus1, ~`depression_binned`)
DESEQ_depression <- DESeq(depression_deseq)
res_depression <- results(DESEQ_depression, tidy=TRUE, 
               #this will ensure that No is your reference group
               contrast = c("depression_binned","Yes","No"))
## Make variable to color by whether it is significant + large change ##
depression_vol_plot <- res_depression %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  ggtitle('Depression DeSeq2 Volcano Plot') +
  theme(plot.title = element_text(hjust = 0.5))+
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))
depression_vol_plot
#### DESeq Volcano Plot for Anxiety ####

anxiety_plus1 <- transform_sample_counts(parkinsons_final_anxiety, function(x) x+1)
anxiety_deseq <- phyloseq_to_deseq2(anxiety_plus1, ~`anxiety_binned`)
DESEQ_anxiety <- DESeq(anxiety_deseq)
res_anxiety <- results(DESEQ_anxiety, tidy=TRUE, 
                          #this will ensure that No is your reference group
                          contrast = c("anxiety_binned","Yes","No"))
## Make variable to color by whether it is significant + large change ##
anxiety_vol_plot <- res_anxiety %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  ggtitle('Anxiety DeSeq2 Volcano Plot') +
  theme(plot.title = element_text(hjust = 0.5))+
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))

#### DESeq Volcano Plot for Anxiety ####

sleep_plus1 <- transform_sample_counts(parkinsons_final_sleep, function(x) x+1)
sleep_deseq <- phyloseq_to_deseq2(sleep_plus1, ~`Sleep_problems`)
DESEQ_sleep <- DESeq(sleep_deseq)
res_sleep <- results(DESEQ_sleep, tidy=TRUE, 
                       #this will ensure that No is your reference group
                       contrast = c("Sleep_problems","Yes","No"))
## Make variable to color by whether it is significant + large change ##
sleep_vol_plot <- res_sleep %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  ggtitle('Sleep Problem DeSeq2 Volcano Plot') +
  theme(plot.title = element_text(hjust = 0.5))+
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))

#### Annotate figures with captions “A, B, C, D” ####


combined_volcano <- plot_grid(depression_vol_plot + theme(legend.position="none"), 
                             anxiety_vol_plot + theme(legend.position="none"), 
                             sleep_vol_plot + theme(legend.position="none"),
                             labels = c('A', 'B','C'), nrow = 1)
legend_for_combined_volcano <- get_legend(
  # create some space to the left of the legend
  depression_vol_plot + theme(legend.box.margin = margin(0, 0, 0, 12))
)

combined_volcano_w_legend <-plot_grid(combined_volcano, legend_for_combined_volcano, rel_widths = c(3, .4))

combined_volcano_w_legend

ggsave('combined_volcano_with_legend.png',
       width = 20,
       height = 10,
       combined_volcano_w_legend)


#!/usr/bin/env Rscript
library(tidyverse)
library(phyloseq)
library(DESeq2)


#### Load data ####
load("~/parkinsons_project_2/DeSeq_PD_Sleep/parkinsons_final_sleep.RData")

head(sam_data(parkinsons_final_sleep))

#### DESeq ####

sp_plus1 <- transform_sample_counts(parkinsons_final_sleep, function(x) x+1)
sp_deseq <- phyloseq_to_deseq2(sp_plus1, ~`Sleep_problems`)
DESEQ_sp <- DESeq(sp_deseq)
res <- results(DESEQ_sp, tidy=TRUE, 
               #this will ensure that No is your reference group
               contrast = c("Sleep_problems","Yes","No"))
View(res)

## Volcano plot: effect size VS significance
sp_vol <- ggplot(res) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj)))
ggsave('sp_deseq.png',
       sp_vol,
       height = 10, width = 10)

## Make variable to color by whether it is significant + large change ##
sp_vol_plot <- res %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))

ggsave(filename="sleep_volcano_plot.png",sp_vol_plot)

# To get table of results
sigASVs <- res %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(sigASVs)
# Get only asv names
sigASVs_vec <- sigASVs %>%
  pull(ASV)


# Prune phyloseq file
sp_DESeq <- prune_taxa(sigASVs_vec,parkinsons_final_sleep)
sigASVs <- tax_table(sp_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

sleep_bar_plot <- ggplot(sigASVs) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
ggsave('sleep_bar_plot.png',
       sleep_bar_plot)

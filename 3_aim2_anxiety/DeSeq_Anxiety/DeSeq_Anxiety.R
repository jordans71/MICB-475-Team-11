#!/usr/bin/env Rscript
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

pkgs <- c("phyloseq", "ALDEx2", "SummarizedExperiment", "Biobase", "devtools", 
          "ComplexHeatmap", "BiocGenerics", "BiocManager", "metagenomeSeq", 
          "Maaslin2", "edgeR", "lefser", "limma", "KEGGREST", "DESeq2")
for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg)
}

library(tidyverse)
library(phyloseq)
library(DESeq2)
library(ggplot2)


#### Load data ####

#### DESeq ####

anxiety_plus1 <- transform_sample_counts(parkinsons_final_anxiety, function(x) x+1)
anxiety_deseq <- phyloseq_to_deseq2(anxiety_plus1, ~`anxiety_binned`)
DESEQ_anxiety <- DESeq(anxiety_deseq)
res <- results(DESEQ_anxiety, tidy=TRUE, 
               #this will ensure that No is your reference group
               contrast = c("anxiety_binned","High","Low"))
View(res)

## Volcano plot: effect size VS significance
anxiety_vol <- ggplot(res) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj)))
ggsave('anxiety_deseq.png',
       anxiety_vol,
       height = 10, width = 10)

## Make variable to color by whether it is significant + large change ##
anxiety_vol_plot <- res %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))
ggsave(filename="anxiety_volcano_plot.png",anxiety_vol_plot)



# To get table of results
sigASVs <- res %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(sigASVs)
# Get only asv names
sigASVs_vec <- sigASVs %>%
  pull(ASV)


# Prune phyloseq file
anxiety_DESeq <- prune_taxa(sigASVs_vec,parkinsons_final_anxiety)
sigASVs <- tax_table(anxiety_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

anxiety_bar_plot <- ggplot(sigASVs) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
anxiety_bar_plot

ggsave('anxiety_bar_plot.png',
       anxiety_bar_plot)

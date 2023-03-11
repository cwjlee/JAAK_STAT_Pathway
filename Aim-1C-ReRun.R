#!/usr/bin/env Rscript

#### Loading packages #### 
library(tidyverse)
library(phyloseq)
library(DESeq2)

#### Loading data #### 
load("/Users/ashleenkhatra/Desktop/MICB_475/Project 2/Re-run/dysautonomia_final.RData")

#### Filtering samples ####
dysautonomia_mild <- subset_samples(dysautonomia_final, Test.FD.severity %in% c("control", "mild"))
dysautonomia_severe <- subset_samples(dysautonomia_final, Test.FD.severity %in% c("control", "severe"))
dysautonomia_mil_severe <- subset_samples(dysautonomia_final, Test.FD.severity %in% c("mild", "severe"))

#### DESeq for Control vs Mild FD ####
deseq_dysautonomia_mild <- phyloseq_to_deseq2(dysautonomia_mild,~Test.FD.severity)
DESEQ_mild <- DESeq(deseq_dysautonomia_mild)
res <- results(DESEQ_mild, tidy=TRUE)
View(res)
# Volcano plot: effect size VS significance
ggplot(res) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj)))
res %>%
  mutate(significant = padj<0.05 & abs(log2FoldChange)>1.5) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))

# Table of results
sigASVs <- res %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(sigASVs)
# Significant ASVs
sigASVs_vec <- sigASVs %>%
  pull(ASV)
# The only significantly different ASV between the control and mild group is 95d6a2a0e7a344010b5fd670677848e1 

# Prune phyloseq file
mild_FD_DESeq <- prune_taxa(sigASVs_vec,dysautonomia_mild)
# Phlyum level comparison
sigASVs <- tax_table(mild_FD_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Phylum = make.unique(Phylum)) %>%
  mutate(Phylum = factor(Phylum, levels=unique(Phylum)))
ggplot(sigASVs) +
  geom_bar(aes(x=Phylum, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Phylum, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE))

# Genus level comparison
sigASVs <- tax_table(mild_FD_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

ggplot(sigASVs) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) 

# Species level comparison
sigASVs <- tax_table(mild_FD_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Species = make.unique(Species)) %>%
  mutate(Species = factor(Species, levels=unique(Species)))

ggplot(sigASVs) +
  geom_bar(aes(x=Species, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Species, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE))

#### DESeq for Control vs Severe FD ####
deseq_dysautonomia_mild_severe <- phyloseq_to_deseq2(dysautonomia_mil_severe,~Test.FD.severity)
DESEQ_severe <- DESeq(deseq_dysautonomia_mild_severe)
res_2 <- results(DESEQ_severe, tidy=TRUE)
View(res_3)
# Volcano plot: effect size VS significance
ggplot(res_3) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj)))
res_3 %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))

# Table of results
sigASVs_severe <- res_3 %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(sigASVs_severe)
# Significant ASVs
sigASVs_severe_vec <- sigASVs_severe %>%
  pull(ASV)
# There are 4 significantly different ASVs between the control and severe group

# Prune phyloseq file
severe_FD_DESeq <- prune_taxa(sigASVs_severe_vec,dysautonomia_severe)
# Phlyum level comparison
severe_sigASVs <- tax_table(severe_FD_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_severe) %>%
  arrange(log2FoldChange) %>%
  mutate(Phylum = make.unique(Phylum)) %>%
  mutate(Phylum = factor(Phylum, levels=unique(Phylum)))
ggplot(severe_sigASVs) +
  geom_bar(aes(x=Phylum, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Phylum, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE))

# Genus level comparison
severe_sigASVs <- tax_table(severe_FD_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_severe) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))
ggplot(severe_sigASVs) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE))

# Species level comparison
severe_sigASVs <- tax_table(severe_FD_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_severe) %>%
  arrange(log2FoldChange) %>%
  mutate(Species = make.unique(Species)) %>%
  mutate(Species = factor(Species, levels=unique(Species)))
ggplot(severe_sigASVs) +
  geom_bar(aes(x=Species, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Species, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE))

#### DESeq for Mild vs Severe FD ####
deseq_dysautonomia_mil_severe <- phyloseq_to_deseq2(dysautonomia_moderate,~FD.severity)
DESEQ_moderate <- DESeq(deseq_dysautonomia_moderate)
res_2 <- results(DESEQ_moderate, tidy=TRUE)
View(res_2)

# Volcano plot: effect size VS significance
ggplot(res_2) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj)))
res_2 %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))

# Table of results
sigASVs_moderate <- res_2 %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
View(sigASVs_moderate)
# There are no ASVs significantly different between the moderate and control groups


#!/usr/bin/env Rscript
library(phyloseq) 
library(ape)
library(tidyverse)
library(vegan)
library(FSA)

#### Loading data #### 
load("/Users/ashleenkhatra/Desktop/MICB_475/Project 2/Re-run2/dysautonomia_final.RData")

#### DESeq for Mild vs Severe FD ####
deseq_dysautonomia_mil_severe <- phyloseq_to_deseq2(dysautonomia_mil_severe,~FD.severity)
DESEQ_mild_severe <- DESeq(deseq_dysautonomia_mil_severe)
res_4 <- results(DESEQ_mild_severe, tidy=TRUE)
View(res_4)

# Volcano plot: effect size VS significance
ggplot(res_4) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj)))
gg_volcano <- res_4 %>%
  mutate(significant = padj<0.05 & abs(log2FoldChange)>1.5) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))

# Table of results
sigASVs_mild_severe <- res_4 %>% 
  filter(padj<0.05 & abs(log2FoldChange)>1.5) %>%
  dplyr::rename(ASV=row)
View(sigASVs_mild_severe)
# Significant ASVs
sigASVs_msevere_vec <- sigASVs_mild_severe %>%
  pull(ASV)
# There are 5 ASVs significantly different between the mild and severe groups: 721bbde09abf2f51fc7d8caeab5b22f1, 0d29a1eb24f264423251f3b2a92f7914, 273fa0191072af3d33e32271b35c8f18, acf76b1f22c9536ca982df6f9b0219da, cdb79eecfa38ac58e99bfdef45fabfce    

# Prune phyloseq file
msevere_FD_DESeq <- prune_taxa(sigASVs_msevere_vec,dysautonomia_mil_severe)
# Phlyum level comparison
msevere_sigASVs <- tax_table(msevere_FD_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_mild_severe) %>%
  arrange(log2FoldChange) %>%
  mutate(Phylum = make.unique(Phylum)) %>%
  mutate(Phylum = factor(Phylum, levels=unique(Phylum)))
ggplot(msevere_sigASVs) +
  geom_bar(aes(x=Phylum, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Phylum, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE))

# Genus level comparison
msevere_sigASVs <- tax_table(msevere_FD_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_mild_severe) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))
ggplot(msevere_sigASVs) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE))+ theme(axis.text.x = element_text(angle = 90))

# Species level comparison
msevere_sigASVs <- tax_table(msevere_FD_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_mild_severe) %>%
  arrange(log2FoldChange) %>%
  mutate(Species = make.unique(Species)) %>%
  mutate(Species = factor(Species, levels=unique(Species)))
ggplot(msevere_sigASVs) +
  geom_bar(aes(x=Species, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Species, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE))+ theme(axis.text.x = element_text(angle = 90))

#### Genus and phyla level comparison between mild and severe conditions ####
# Combining p-values 
combine_pvalues = function(p){
  return(1-pchisq(-2*sum(log(p),na.rm=T),2*sum(!is.na(p))))}
#Running DESEQ on the phyloseq object
deseq_dysautonomia_mil_severe <- phyloseq_to_deseq2(dysautonomia_mil_severe,~FD.severity)
DESEQ_mild_severe <- DESeq(deseq_dysautonomia_mil_severe)
res_4 <- results(DESEQ_mild_severe, tidy=TRUE)
#Combining the results output with taxa information.
tax_ms_dysautoomia <- tax_table(dysautonomia_mil_severe)
dysautonomia_DEseq_results <- data.frame(cbind(res_4,tax_ms_dysautoomia)) 
# Key step is adding the taxa info to the deseq results. Need to make sure the number of rows is the same.
# Filtering the deseq results at different taxa levels.
dysautonomia_DEseq_results_genus_summary<-dysautonomia_DEseq_results %>%
  group_by("Phylum","Class","Order","Family","Genus") %>%
  select("Phylum","Class","Order","Family","Genus",log2FoldChange,padj,baseMean,pvalue) %>%
  dplyr::summarise(mean_2_fold_change=mean(log2FoldChange),padj=combine_pvalues(padj),countmean = mean(baseMean),pvalue=combine_pvalues(pvalue))

# Filtering the deseq results at different the family level
DEseq_results_family_summary<-dysautonomia_DEseq_results %>%
  group_by("Class","Order","Family","Genus") %>%
  select("Class","Order","Family","Genus",log2FoldChange,padj,baseMean,pvalue) %>%
  dplyr::summarise(mean_2_fold_change=mean(log2FoldChange),padj=combine_pvalues(padj),countmean = mean(baseMean),pvalue=combine_pvalues(pvalue))

# Filtering the deseq results at different the order level
DEseq_results_order_summary<-dysautonomia_DEseq_results %>%
  group_by("Order","Family","Genus") %>%
  select("Order","Family","Genus",log2FoldChange,padj,baseMean,pvalue) %>%
  dplyr::summarise(mean_2_fold_change=mean(log2FoldChange),padj=combine_pvalues(padj),countmean = mean(baseMean),pvalue=combine_pvalues(pvalue))

# Filtering the deseq results at different the family level
DEseq_results_class_summary<-dysautonomia_DEseq_results %>%
  group_by("Family","Genus") %>%
  select("Family","Genus",log2FoldChange,padj,baseMean,pvalue) %>%
  dplyr::summarise(mean_2_fold_change=mean(log2FoldChange),padj=combine_pvalues(padj),countmean = mean(baseMean),pvalue=combine_pvalues(pvalue))

# Filtering the deseq results at different the genus level
DEseqafterfiltergenus <- ggplot(dysautonomia_DEseq_results_genus_summary,aes(x=reorder("Phylum", -mean_2_fold_change),y=mean_2_fold_change,color = "Genus"))+
  geom_point() +
  #geom_hline(yintercept = 0,aes(linewidth=1))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = -90,hjust = 0,vjust = 0.5),
        text = element_text(face = "bold"),
        axis.text.x.bottom  = element_text(,angle = -45),)+
  scale_color_manual(values =  c("#2F3CBE",                                       
                                          "#0dc34e",                                                                   
                                          "#e24223",
                                          "#86007D"))+
                                            labs(x="Genus",y="Log2 Fold Change",color= "Phylum")
+ ggtitle("DESeq")


write.table(PAE_DEseq_results_genus_summary,file = "Genus list for females.txt")
ggsave(paste(dataset_name,"-DEseq after filtering PAEvs control_female (genus).png"),plot= DEseqafterfiltergenus)

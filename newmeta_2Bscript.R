#load packages
library(dplyr)
library(phyloseq)
library(ape)
library(tidyverse)
library(vegan)

#load data
metafp <- "new_metadata.tsv"
meta <- read_delim(metafp, delim="\t")

otufp <- "feature-table.txt"
otu <- read_delim(file = otufp, delim="\t", skip=1)

taxfp <- "taxonomy.tsv"
tax <- read_delim(taxfp, delim="\t")

phylotreefp <- "tree.nwk"
phylotree <- read.tree(phylotreefp)

#generate table of disease severity by sex
meta_filter <- select(meta, Sex, Test.FD.severity)
counts <- table(meta_filter$Sex, meta_filter$Test.FD.severity)
counts

#format data for phyloseq object

samp_df <- as.data.frame(meta[,-1]) 
rownames(samp_df)<- meta$`sample-id`
SAMP <- sample_data(samp_df)

otu_mat <- as.matrix(otu[,-1]) 
rownames(otu_mat) <- otu$`#OTU ID`
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) 

tax_mat <- tax %>% 
  select(-Confidence) %>%
  separate(col=Taxon, sep="; ", into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix()

rownames(tax_mat) <- tax$`Feature ID`
TAX <- tax_table(tax_mat)

dysphylo <- phyloseq(OTU, SAMP, TAX, phylotree)
dysphylo_final <- na.omit(dysphylo)
save(dysphylo_final, file = "dysphylo.RData")

#plot beta diversity metrics
bdm <- distance(dysphylo_final, method = "bray")
jdm <- distance(dysphylo_final, method = "jaccard", binary = TRUE)
wufdm <- distance(dysphylo_final, method = "wunifrac")
ufdm <- distance(dysphylo_final, method = "unifrac")

BC <- ordinate(physeq = dysphylo_final, method = "PCoA", distance = bdm)
JCC <- ordinate(physeq = dysphylo_final, method = "PCoA", distance = jdm)
WUF <- ordinate(physeq = dysphylo_final, method = "PCoA", distance = wufdm)
UF <- ordinate(physeq = dysphylo_final, method = "PCoA", distance = ufdm)

BC_plot <- plot_ordination(dysphylo_final, BC, color = "Test.FD.severity", shape = "Sex", title = "Bray-Curtis")
JCC_plot <- plot_ordination(dysphylo_final, JCC, color = "Test.FD.severity", shape = "Sex", title = "Jaccard")
WUF_plot <- plot_ordination(dysphylo_final, WUF, color = "Test.FD.severity", shape = "Sex", title = "Weighted Unifrac")
UF_plot <- plot_ordination(dysphylo_final, UF, color = "Test.FD.severity", shape = "Sex", title = "Unifrac")

print(BC_plot)
print(JCC_plot)
print(WUF_plot)
print(UF_plot)

#statistical testing with PERMANOVA

meta_filter_na <- na.omit(meta_filter)
adonis2(bdm ~ Test.FD.severity*Sex, data = meta_filter_na)
adonis2(jdm ~ Test.FD.severity*Sex, data = meta_filter_na)
adonis2(wufdm ~ Test.FD.severity*Sex, data = meta_filter_na)
adonis2(ufdm ~ Test.FD.severity*Sex, data = meta_filter_na)



#!/usr/bin/env Rscript

#### Loading packages #### 
library(phyloseq) 
library(ape)
library(tidyverse)
library(vegan)
library(FSA)

#### Loading data #### 
metafp <- "/Users/ashleenkhatra/Desktop/MICB_475/Project 2/Re-run/new_metadata.tsv"
meta <- read_delim(metafp, delim="\t")

otufp <- "/Users/ashleenkhatra/Desktop/MICB_475/Project 2/Re-run/feature-table.txt"
otu <- read_delim(file = otufp, delim="\t", skip=1)

taxfp <- "/Users/ashleenkhatra/Desktop/MICB_475/Project 2/Re-run/taxonomy.tsv"
tax <- read_delim(taxfp, delim="\t")

phylotreefp <- "/Users/ashleenkhatra/Desktop/MICB_475/Project 2/Re-run/tree.nwk"
phylotree <- read.tree(phylotreefp)

#### Formating the OTU table for phyloseq ####
# saving everything except first column (OTU ID) into a matrix
otu_mat <- as.matrix(otu[,-1])
# Making the first column (#OTU ID) the rownames of the new matrix
rownames(otu_mat) <- otu$'#OTU ID'
# Using the "otu_table" function to make an OTU table
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) 

#### Formating the sample metadata table for phyloseq ####
# Saving everything except sampleid as new data frame
samp_df <- as.data.frame(meta[,-1])
# Making the sampleids the rownames
rownames(samp_df)<- meta$"sample-id"
# Making phyloseq sample data with sample_data() function
SAMP <- sample_data(samp_df)

#### Formating the taxonomy table ####
# Converting the taxon strings to a table with separate taxa rank columns
tax_mat <- tax %>% select(-Confidence)%>%
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix()
# Saving everything except feature IDs 
tax_mat <- tax_mat[,-1]
# Making the sampleids the rownames
rownames(tax_mat) <- tax$`Feature ID`
# Making a taxa table
TAX <- tax_table(tax_mat)

#### Create phyloseq object ####
# Merge all into a phyloseq object
dysautonomia <- phyloseq(OTU, SAMP, TAX, phylotree)

#### Looking at phyloseq object #####
# Viewing components of phyloseq object
otu_table(dysautonomia)
sample_data(dysautonomia)
tax_table(dysautonomia)
phy_tree(dysautonomia)

#### Filtering the data ####
# Removing chloroplasts, mitochondria, non-bacterial samples
dysautonomia_filt <- subset_taxa(dysautonomia,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")
# Removing the ASVs that have less than 5 counts total (low abundance ASVs)
dysautonomia_filt_nolow <- filter_taxa(dysautonomia_filt, function(x) sum(x)>5, prune = TRUE)
# Removing low quality samples with less than 100 reads
dysautonomia_filt_nolow_samps <- prune_samples(sample_sums(dysautonomia_filt_nolow)>100, dysautonomia_filt_nolow)
# Remove samples where Test.FD.severity is na
dysautonomia_final <- subset_samples(dysautonomia_filt_nolow_samps, !is.na(Test.FD.severity) )


#### Rarefy samples ####
rarecurve(t(as.data.frame(otu_table(dysautonomia_final))), cex=0.1)
dysautonomia_rare <- rarefy_even_depth(dysautonomia_final, rngseed = 1, sample.size = 9000)
# 54 out of 214 samples were removed
# 80 out of ____ OTUs were removed

#### Saving the rarefy samples ####
save(dysautonomia_final, file="/Users/ashleenkhatra/Desktop/MICB_475/Project 2/Re-run/dysautonomia_final.RData")
save(dysautonomia_rare, file="/Users/ashleenkhatra/Desktop/MICB_475/Project 2/Re-run/dysautonomia_rare.RData")

#### Load in RData ####
load("/Users/ashleenkhatra/Desktop/MICB_475/Project 2/Re-run/dysautonomia_final.RData")
load("/Users/ashleenkhatra/Desktop/MICB_475/Project 2/Re-run/dysautonomia_rare.RData")

#### Alpha diversity metrics ####
# Shannon and Chao1 diversity plots
gg_richness <- plot_richness(dysautonomia_rare, x = "Test.FD.severity", measures = c("Shannon","Chao1")) +
  xlab("FD Severity Level") +
  geom_boxplot()

gg_shannon <- plot_richness(dysautonomia_rare, x = "Test.FD.severity", measures = c("Shannon")) +
  xlab("FD Severity Level") +
  geom_boxplot()

gg_chao1 <- plot_richness(dysautonomia_rare, x = "Test.FD.severity", measures = c("Chao1")) +
  xlab("FD Severity Level") +
  geom_boxplot()

# Statistical analysis for alpha diversity metrics
# Shannon diversity Kruskall-Wallis Test
alphadiv <- estimate_richness(dysautonomia_rare)
samp_dat <- sample_data(dysautonomia_rare)
samp_dat_wdiv <- data.frame(samp_dat, alphadiv)

ggplot(samp_dat_wdiv) +
  geom_point(aes(x=Test.FD.severity, y=Shannon))

kruskal.test(Shannon ~ Test.FD.severity, data = samp_dat_wdiv)
?kruskal.test
# p-value = 0.577
# Dunn test for pairwise comparison
dunnTest(Shannon ~ Test.FD.severity, data = samp_dat_wdiv, method="bh")   

# Chao1 diversity Kruskall-Wallis Test
kruskal.test(Chao1 ~ Test.FD.severity, data = samp_dat_wdiv)
# p-value = 0.001743
# Dunn test for pairwise comparison
dunnTest(Chao1 ~ Test.FD.severity, data = samp_dat_wdiv, method="bh")   

# Saving the alpha diversity bar plots
ggsave(filename = "/Users/ashleenkhatra/Desktop/MICB_475/Project 2/Re-run/alpha_diversity_richness.png"
       , gg_richness
       , height=4, width=6)

ggsave(filename = "/Users/ashleenkhatra/Desktop/MICB_475/Project 2/Re-run/shannon_diversity_richness.png"
       , gg_shannon
       , height=4, width=6)

ggsave(filename = "/Users/ashleenkhatra/Desktop/MICB_475/Project 2/Re-run/Chao1_diversity_richness.png"
       , gg_chao1
       , height=4, width=6)


#### Beta diversity metrics ####
# Bray Curtis
bc_dm <- distance(dysautonomia_rare, method="bray")
pcoa_bc <- ordinate(dysautonomia_rare, method="PCoA", distance=bc_dm)
plot_ordination(dysautonomia_rare, pcoa_bc, color = "Test.FD.severity") + labs(col = "FD Severity Level") + stat_ellipse()
# Bray Curtic diversity PERMANOVA Test
dm_bray <- vegdist(t(otu_table(dysautonomia_rare)), method="bray")
adonis2(dm_bray ~ Test.FD.severity, data=samp_dat_wdiv)
# p-value = 0.001

# Weighted Unifrac
bc_dm <- distance(dysautonomia_rare, method="wunifrac")
pcoa_bc <- ordinate(dysautonomia_rare, method="PCoA", distance=bc_dm)
plot_ordination(dysautonomia_rare, pcoa_bc, color = "Test.FD.severity") + labs(col = "FD Severity Level") + stat_ellipse()
# Weighted Unifrac Curtic diversity PERMANOVA Test
dm_unifrac <- UniFrac(dysautonomia_rare, weighted=TRUE)
adonis2(dm_unifrac ~ Test.FD.severity, data=samp_dat_wdiv)
# p-value = 0.001

# Jaccard
bc_dm <- distance(dysautonomia_rare, method = "jaccard", binary = TRUE)
pcoa_bc <- ordinate(dysautonomia_rare, method="PCoA", distance=bc_dm)
plot_ordination(dysautonomia_rare, pcoa_bc, color = "Test.FD.severity") + labs(col = "FD Severity Level") + stat_ellipse()
# Jaccard diversity PERMANOVA Test
dm_jaccard <- vegdist(t(otu_table(dysautonomia_rare)), method="jaccard")
adonis2(dm_jaccard ~ Test.FD.severity, data=samp_dat_wdiv)

# Unweighted Unifrac
bc_dm <- distance(dysautonomia_rare, method="unifrac")
pcoa_bc <- ordinate(dysautonomia_rare, method="PCoA", distance=bc_dm)
plot_ordination(dysautonomia_rare, pcoa_bc, color = "Test.FD.severity")
gg_unifrac <- plot_ordination(dysautonomia_rare, pcoa_bc, color = "Test.FD.severity") +
  labs(col = "FD Severity Level") + stat_ellipse()

ggsave(filename = "Unweighted_Unifrac_Beta_Diversity.png"
       , gg_unifrac
       , height=4, width=6)
# (Unweighted) Unifrac diversity PERMANOVA Test
dm_unifrac <- UniFrac(dysautonomia_rare, weighted=FALSE)
adonis2(dm_unifrac ~ Test.FD.severity, data=samp_dat_wdiv)
# p-value = 0.001


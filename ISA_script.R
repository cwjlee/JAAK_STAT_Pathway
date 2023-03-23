

library(tidyverse)
library(phyloseq)
library(indicspecies)
library(ape) # importing trees
library(vegan)

#### Load data ####
# Change file paths as necessary
metafp <- "/Users/jennyshee/Desktop/W23/MICB475/ISA/new_metadata.tsv"
meta <- read_delim(metafp, delim="\t")

otufp <- "/Users/jennyshee/Desktop/W23/MICB475/JAAK_STAT_project/dysautonomia/exported_data/Output/table_rmmcna_export/feature-table.txt"
otu <- read_delim(file = otufp, delim="\t", skip=1)

taxfp <- "/Users/jennyshee/Desktop/W23/MICB475/JAAK_STAT_project/dysautonomia/exported_data/Output/taxonomy_export/taxonomy.tsv"
tax <- read_delim(taxfp, delim="\t")

phylotreefp <- "/Users/jennyshee/Desktop/W23/MICB475/JAAK_STAT_project/dysautonomia/exported_data/Output/rooted_tree_export/tree.nwk"
phylotree <- read.tree(phylotreefp)

#### Format OTU table ####
# save everything except first column (OTU ID) into a matrix
otu_mat <- as.matrix(otu[,-1])
# Make first column (#OTU ID) the rownames of the new matrix
rownames(otu_mat) <- otu$`#OTU ID`
# Use the "otu_table" function to make an OTU table
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) 

#### Format sample metadata ####
# Save everything except sampleid as new data frame
samp_df <- as.data.frame(meta[,-1])
# Make sampleids the rownames
rownames(samp_df)<- meta$sampleid
# Make phyloseq sample data with sample_data() function
SAMP <- sample_data(samp_df)

#### Formatting taxonomy ####
# Convert taxon strings to a table with separate taxa rank columns
tax_mat <- tax %>% 
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix() # Saving as a matrix
# Save everything except feature IDs 
tax_mat <- tax_mat[,-1]
# Make sampleids the rownames
rownames(tax_mat) <- tax$`Feature ID`
# Make taxa table
TAX <- tax_table(tax_mat)

#### Create phyloseq object ####
# Merge all into a phyloseq object
meta_ps <- phyloseq(OTU, SAMP, TAX, phylotree)

#### ISA Phylum ####
meta_phylum <- tax_glom(meta_ps, "Phylum", NArm = FALSE)
meta_phylum_RA <- transform_sample_counts(meta_phylum, fun=function(x) x/sum(x))
#ISA
isa_meta_phylum <- multipatt(t(otu_table(meta_phylum_RA)), cluster = sample_data(meta_phylum_RA)$FD.severity)
summary(isa_meta_phylum)
taxtable <- tax_table(meta_ps) %>% as.data.frame() %>% rownames_to_column(var="ASV")

isa_phylum <- isa_meta_phylum$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable) %>%
  filter(p.value<0.05) %>% as.data.frame()

#### ISA Class ####
meta_class <- tax_glom(meta_ps, "Class", NArm = FALSE)
meta_class_RA <- transform_sample_counts(meta_class, fun=function(x) x/sum(x))
#ISA
isa_meta_class <- multipatt(t(otu_table(meta_class_RA)), cluster = sample_data(meta_class_RA)$FD.severity)
summary(isa_meta_class)
taxtable <- tax_table(meta_ps) %>% as.data.frame() %>% rownames_to_column(var="ASV")


isa_class <- isa_meta_class$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable) %>%
  filter(p.value<0.05) %>% as.data.frame()

#### ISA Order ####
meta_order <- tax_glom(meta_ps, "Order", NArm = FALSE)
meta_order_RA <- transform_sample_counts(meta_order, fun=function(x) x/sum(x))
#ISA
isa_meta_order <- multipatt(t(otu_table(meta_order_RA)), cluster = sample_data(meta_order_RA)$FD.severity)
summary(isa_meta_order)
taxtable <- tax_table(meta_ps) %>% as.data.frame() %>% rownames_to_column(var="ASV")


isa_order <- isa_meta_order$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable) %>%
  filter(p.value<0.05) %>% as.data.frame()

#### ISA Family ####
meta_family <- tax_glom(meta_ps, "Family", NArm = FALSE)
meta_family_RA <- transform_sample_counts(meta_family, fun=function(x) x/sum(x))
#ISA
isa_meta_family <- multipatt(t(otu_table(meta_family_RA)), cluster = sample_data(meta_family_RA)$FD.severity)
summary(isa_meta_family)
taxtable <- tax_table(meta_ps) %>% as.data.frame() %>% rownames_to_column(var="ASV")


isa_family <- isa_meta_family$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable) %>%
  filter(p.value<0.05) %>% as.data.frame()

#### ISA Genus ####
meta_genus <- tax_glom(meta_ps, "Genus", NArm = FALSE)
meta_genus_RA <- transform_sample_counts(meta_genus, fun=function(x) x/sum(x))
#ISA
isa_meta_genus <- multipatt(t(otu_table(meta_genus_RA)), cluster = sample_data(meta_genus_RA)$FD.severity)
summary(isa_meta_genus)
taxtable <- tax_table(meta_ps) %>% as.data.frame() %>% rownames_to_column(var="ASV")

isa_genus <- isa_meta_genus$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable) %>%
  filter(p.value<0.05) %>% as.data.frame()

#### ISA Species ####
meta_species <- tax_glom(meta_ps, "Species", NArm = FALSE)
meta_species_RA <- transform_sample_counts(meta_species, fun=function(x) x/sum(x))
#ISA
isa_meta_species <- multipatt(t(otu_table(meta_species_RA)), cluster = sample_data(meta_species_RA)$FD.severity)
summary(isa_meta_genus)
taxtable <- tax_table(meta_ps) %>% as.data.frame() %>% rownames_to_column(var="ASV")

isa_species <- isa_meta_species$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable) %>%
  filter(p.value<0.05) %>% as.data.frame() 
  

# visualize


  

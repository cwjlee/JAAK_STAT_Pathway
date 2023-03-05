##### 2B-1 load packages #####
library(phyloseq)
library(ape)
library(tidyverse)
library(vegan)

##### 2B-2 load exported data files (sample-data.txt, feature-table.txt, taxonomy.tsv, tree.nwk) #####
metafp <- "/Users/jennyshee/Desktop/W23/MICB475/JAAK_STAT_project/dysautonomia/metadata.tsv"
meta <- read_delim(metafp, delim="\t")

otufp <- "/Users/jennyshee/Desktop/W23/MICB475/JAAK_STAT_project/dysautonomia/exported_data/Output/table_rmmcna_export/feature-table.txt"
otu <- read_delim(file = otufp, delim="\t", skip=1)

taxfp <- "/Users/jennyshee/Desktop/W23/MICB475/JAAK_STAT_project/dysautonomia/exported_data/Output/taxonomy_export/taxonomy.tsv"
tax <- read_delim(taxfp, delim="\t")

phylotreefp <- "/Users/jennyshee/Desktop/W23/MICB475/JAAK_STAT_project/dysautonomia/exported_data/Output/rooted_tree_export/tree.nwk"
phylotree <- read.tree(phylotreefp)

##### 2B-3 View sample metadata and generate histogram of pathology severity by sex groups of male and female to check for sufficient sample to conduct subsequent analysis #####
meta_FD.severity_Sex <- select(meta, FD.severity, Sex)
graph_without_count <- ggplot(meta_FD.severity_Sex, aes(x=FD.severity, y=Sex)) + geom_count() 
graph_with_count <- graph_without_count + geom_text(data = ggplot_build(graph_without_count)$data[[1]],aes(x, y, label = n, vjust= -2))

##### 2B-4 Format sample metadata, OTU table, and taxonomy. Combine and save all as a phyloseq object (output = dysautonomia_final.RData) #####
#### Format metadata ####
# Save everything except sampleid as new data frame
samp_df <- as.data.frame(meta[,-1]) #remove first column, which are samples IDs
# Make sampleids the rownames
rownames(samp_df)<- meta$sampleid
# Make phyloseq sample data with sample_data() function
SAMP <- sample_data(samp_df)

#### Format OTU table ####
# save everything except first column (OTU ID) into a matrix
otu_mat <- as.matrix(otu[,-1]) #don't want first column, so remove by otu[,1]
# Make first column (#OTU ID) the rownames of the new matrix
rownames(otu_mat) <- otu$`#OTU ID`
# Use the "otu_table" function to make an OTU table
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) 

#### Format taxonomy ####
# Convert taxon strings to a table with separate taxa rank columns
tax_mat <- tax %>% select(-Confidence)%>%
  separate(col=Taxon, sep="; " 
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix() 
# Save everything except feature IDs 
tax_mat <- tax_mat[,-1]
# Make sampleids the rownames
rownames(tax_mat) <- tax$`Feature ID`
# Make taxa table
TAX <- tax_table(tax_mat)

#### Create phyloseq object ####
# Merge all into a phyloseq object
dysautonomia_final <- phyloseq(OTU, SAMP, TAX, phylotree)
# Save phylosec object
save(dysautonomia_final, file="dysautonomia_final.RData")

##### 2B-5 Plot beta diversity using distance, ordinate and plot_ordination commands where respective methods are specified. Generate plots for Bray-Curtis, Weighted and unweighted Unifrac, and Jaccard #####
?distance
?ordinate
?plot_ordination
# unweighted unifrac
distance(physeq = dysautonomia_final, method = "unifrac")
# weighted unifrac
distance(physeq = dysautonomia_final, method = "wunifrac")
# jaccard 
distance(physeq = dysautonomia_final, method = "jaccard", binary = TRUE)
# Bray-Curtis
BC <- ordinate(physeq = dysautonomia_final, method = "PCoA", distance = "bray")

BC_plot = plot_ordination(dysautonomia_final, BC, color = "FD.severity", shape = "Sex", title = "Bray-Curtis Plot")
print(BC_plot)


##### 2B-6 Perform statistical analysis for beta diversity analysis using R #####












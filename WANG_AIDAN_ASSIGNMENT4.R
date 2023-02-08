#!/usr/bin/env RScript

#### ASSIGNMENT 4 ####
# For each question, write the code necessary to produce the answer.
# All answers **MUST OUTPUT THE SOLUTION DIRECTLY**

# E.g. if a question asks you how many transects there are,
# you cannot print(meta) and then manually count the number of sites on
# the screen. You must use R functions to count the number of sites for you
# so that the output in the console is the number of sites.
# --> This is not acceptable
print(meta$transect.name) # I have counted two transects
# --> This is acceptable
length(unique(meta$transect.name))

### Load in the atacama soil sample metadata and OTU table

metaFP <- "atacamasoil_sample-data.txt"
meta <- read.delim(file=metaFP)

otuFP <- "atacamasoil_feature-table.txt"
otu <- read.delim(file=otuFP, skip=1, row.names = 1)

#### Q1: How many samples have vegetation and no vegetation ("vegetation") in each transect? (ie, how many samples are in each vegetation-transect pair)  ####

# Extract unique values
vegetation <- unique(meta$vegetation)
transectnames <- unique(meta$transectname)

# Loop over unique values
for (t in transectnames) {
  for (v in vegetation) {
    count <- sum(meta$transectname == t & meta$vegetation == v)
    print(paste("The number of samples belonging to", t, "with a vegetation value of", v, "is", count))
  }
}


#### Q2: Create a loop that will print the transect name, then print the number of rows that belong to each transect  ####
# Hint: You will need to filter the metadata by transect first, then calculate the number of rows in the loop

#filter metadata by transect name
meta_trans <- meta[,7]

# Extract unique values
transects <- unique(meta_trans)

# Loop over unique values
for (u in transects) {
  count <- sum(meta_trans == u)
  print(paste("There are", count, "samples belonging to", u))
}


#### Q3: What sample has the smallest number of reads in the OTU table? ####

#determine minimum read length and sample ID

sort(colSums(otu), decreasing = FALSE)[1]

#Sample ID with the lowest read length is YUN3259.1.1, with a read length of 6



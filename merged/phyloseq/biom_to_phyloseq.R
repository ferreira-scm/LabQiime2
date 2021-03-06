setwd("/SAN/Susanas_den/gitProj/LabQiime2/merged/phyloseq/")

library(biomformat)
library(phyloseq)
library(ggplot2)
library(ape)
library(tidyr)

source("/SAN/Susanas_den/EimeriaMicrobiome/R/1_Data_preparation.R")

otu <- read.table("otu_table.txt", header=TRUE)
tax <- read.table("taxonomy.tsv", sep = "\t", header = TRUE)


rownames(otu) <- otu[,1]
otu <- otu[,-1]
tax[1,1]
rownames(tax) <- tax[,1]
tax <- tax[,-1]
tax[1,]

tax <- tax %>%
    separate(Taxon, c("domain", "phyla", "class", "order", "family", "genus", "species"),"; ")

head(tax)

OTU <- otu_table(as.matrix(otu), taxa_are_rows = TRUE)
TAX <- tax_table(as.matrix(tax))
META <- sample_data(sample.data)

sample_names(OTU)  <- gsub(".*\\.", "", sample_names(OTU))
sample_names(OTU)
sample_names(META)
taxa_names(OTU)
taxa_names(TAX)

ps_d <- phyloseq(OTU, TAX, META)

saveRDS(ps_d, file="/SAN/Susanas_den/gitProj/LabQiime2/merged/tmp/deblur_phyloseq_merged.RDS")

rank_names(ps_d)
sum(otu_table(subset_taxa(ps_d, domain%in%"d__Eukaryota")))/sum(otu_table(ps_d))
sum(otu_table(subset_taxa(ps_d, family%in%"f__Eimeriorina")))/sum(otu_table(ps_d))

#filter by family
ps_d_eimf <- subset_taxa(ps_d, family%in%"f__Eimeriorina")
ps_d_eimg <- subset_taxa(ps_d, genus%in%"g__Eimeria")

setwd("/SAN/Susanas_den/gitProj/LabQiime2/sorted/phyloseq/")

library(biomformat)
library(phyloseq)
library(ggplot2)
library(ape)
source("/SAN/Susanas_den/EimeriaMicrobiome/R/1_Data_preparation.R")
otu <- read.table("otu_table.txt", header=TRUE)
tax <- read.table("taxonomy.tsv", sep = "\t", header = TRUE)

rownames(otu) <- otu[,1]
otu <- otu[,-1]

tax[1,1]
rownames(tax) <- tax[,1]
tax <- tax[,-1]
tax[1,]

library(tidyr)

tax <- tax %>%
    separate(Taxon, c("domain", "phyla", "class", "order", "family", "genus", "species"),"; ")


OTU <- otu_table(as.matrix(otu), taxa_are_rows = TRUE)
TAX <- tax_table(as.matrix(tax))
META <- sample_data(sample.data)
sample_names(OTU)  <- gsub(".*\\.", "", sample_names(OTU))

sample_names(OTU)
sample_names(META)
taxa_names(OTU)
taxa_names(TAX)

ps_ma_d_18S <- phyloseq(OTU, TAX, META)

rank_names(ps_ma_d_18S)
sum(otu_table(subset_taxa(ps_ma_d_18S, domain%in%"d__Eukaryota")))/sum(otu_table(ps_ma_d_18S))
sum(otu_table(subset_taxa(ps_ma_d_18S, family%in%"f__Eimeriorina")))/sum(otu_table(ps_ma_d_18S))

#filter by family
ps_ma_d_18S_eimf <- subset_taxa(ps_d_18S, family%in%"f__Eimeriorina")

saveRDS(ps_ma_d_18S, file="/SAN/Susanas_den/gitProj/LabQiime2/sorted/tmp/18SmaSorted_deblur_phyloseq.RDS")

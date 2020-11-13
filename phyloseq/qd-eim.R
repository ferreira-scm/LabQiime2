library(biomformat)

library(phyloseq)


#feature_table <- system.file("extdata", "feature-table.biom", package="biomformat")

#feature <- read_biom(feature_table)

#feature <- read_biom("feature-table.biom")

otu <- read.table("otu_table.txt", header=TRUE)

head(otu)

tax <- read.table("taxonomy.tsv", sep = "\t", header = TRUE)

head(tax)

merged_file <- merge(otu, tax, by.x=c("OTUID"), by.y=("OTUID"))

head(merged_file)

write.table(merged_file, "combined_otu_tax.tsv", sep="\t", col.names=TRUE, row.names=FALSE)


library(ggplot2)

library(ape)

source("/SAN/Susanas_den/EimeriaMicrobiome/R/1_Data_preparation.R")

rownames(otu) <- otu[,1]
otu <- otu[,-1]

head(tax)


tax[1,1]

rownames(tax) <- tax[,1]

tax <- tax[,-1]

tax[1,]

library(tidyr)

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

ps <- phyloseq(OTU, TAX, META)


rank_names(ps)

get_taxa_unique(ps)

sum(otu_table(subset_taxa(ps, domain%in%"d__Eukaryota")))/sum(otu_table(ps))

sum(otu_table(subset_taxa(ps, family%in%"f__Eimeriorina")))/sum(otu_table(ps))

                                        #filter by family

ps_eimf <- subset_taxa(ps, family%in%"f__Eimeriorina")


ps_eimg <- subset_taxa(ps, genus%in%"g__Eimeria")


sdtmulti <- data.frame(sample_sums(otu_table(ps_eimf)))

sdtmulti[,2] <- rownames(sdtmulti)

head(sdtmulti)

sdtmulti[,3] <- as.data.frame(sample_sums(ps_eimg))

sdtmulti[,4] <- as.data.frame(sample_sums(ps))

colnames(sdtmulti) <- c("Reads_Eimeriorina", "labels", "Reads_Eimeriidae", "Reads_total")



sdt<- read.csv(file = "/SAN/Susanas_den/EimeriaMicrobiome/R/results/sdt.csv")

if(!exists("data.inf.exp")){
  data.inf.exp<- read.csv(file="/SAN/Victors_playground/Eimeria_microbiome/qPCR/sample_data_qPCR.csv")
}


##Keep useful information

keeps <- c("labels", "TotalReads", "ReadsEim")
sdt <- sdt[,colnames(sdt) %in%keeps]

##Get unique labels from qPCR data
keeps <- c("labels", "Qty_mean", "Genome_copies_mean","Tm_mean", "Infection")
data.inf.exp <- data.inf.exp[,colnames(data.inf.exp)%in%keeps]

### Join all the data in the same dataframe 
sdt<- join(sample.data, sdt, by="labels") ## First sample data and multimarker read data 
sdt<- join(sdt, data.inf.exp, by="labels") ## then qPCR data
sdt <- join(sdt, sdtmulti, by="labels")

head(sdt)

sdt$dpi<- as.factor(sdt$dpi)

opg_multi <-  ggplot(sdt, aes(log(1+OPG), log(1+Reads_Eimeriorina)))+
  geom_jitter(shape=21, position=position_jitter(0.2), size=4, aes(fill= dpi), color= "black", alpha=0.7)+
  labs(tag= "B)")+
  theme_bw()+
  theme(text = element_text(size=16))

opg_qpcr <- ggplot(sdt, aes(log(1+OPG), log(1+Genome_copies_mean)))+
  geom_jitter(shape=21, position=position_jitter(0.2), size=4, aes(fill= dpi), color= "black", alpha=0.7)+
  labs(tag= "A)")+
  theme_bw()+
  theme(text = element_text(size=16))

qpcr_multi <-  ggplot(sdt, aes(log(1+Genome_copies_mean), log(1+Reads_Eimeriorina)))+
  geom_jitter(shape=21, position=position_jitter(0.2), size=4, aes(fill= dpi), color= "black", alpha=0.7)+
  labs(tag= "C)")+
  theme_bw()+
  theme(text = element_text(size=16))


opg_qpcr

opg_multi

qpcr_multi


install.packages("gridExtra")

library(gridExtra)

pdf("dna_comp.pdf", width=10, height=7)
grid.arrange(opg_qpcr, opg_multi, qpcr_multi, ncol=2, nrow=2)
dev.off()

str(sdt$dpi)

dpi_opg <- ggplot(sdt, aes(as.numeric(dpi)), log(1+OPG), colour= EH_ID))+
 # xlab("Day post infection")+
 # scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9,10))+
  geom_jitter(shape=21, position=position_jitter(0.0), size=2.5, aes(fill= EH_ID), color= "black", alpha= 0.7)


dpi_multi <- ggplot(sdt, aes(dpi), log(1+Reads_Eimeriorina), colour= EH_ID))+
  xlab("Day post infection")+
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9,10))+
  geom_jitter(shape=21, position=position_jitter(0.0), size=2.5, aes(fill= EH_ID), color= "black", alpha= 0.7)+
  labs(tag= "B)")+
  theme_bw()+
  theme(text = element_text(size=16))

dpi_qpcr <- ggplot(sdt, aes(dpi, log(1+Genome_copies_mean), colour= (EH_ID)))+
                        geom_jitter(shape=21, position=position_jitter(0.0), size=2.5, aes(fill= EH_ID), color= "black", alpha= 0.7) +
                        #xlab("Day post infection")+
                        #scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9,10))+
                        labs(tag= "C)")+
                        theme_bw()+
                        theme(text = element_text(size=16))


dpi_opg

dpi_multi

dpi_qpcr

pdf("dpi_comp.pdf", width=10, height=7)
grid.arrange(dpi_opg, dpi_multi, dpi_qpcr, ncol=2, nrow=2)
dev.off()

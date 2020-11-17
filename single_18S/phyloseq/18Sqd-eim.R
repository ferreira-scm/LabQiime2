setwd("/SAN/Susanas_den/gitProj/LabQiime2/single_18S/phyloseq/")

library(biomformat)
library(phyloseq)
library(ggplot2)
library(ape)
source("/SAN/Susanas_den/EimeriaMicrobiome/R/1_Data_preparation.R")
otu <- read.table("otu_table.txt", header=TRUE)

head(otu)

tax <- read.table("taxonomy.tsv", sep = "\t", header = TRUE)

head(tax)

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
ps <- phyloseq(OTU, TAX, META)


rank_names(ps)
get_taxa_unique(ps)
sum(otu_table(subset_taxa(ps, domain%in%"d__Eukaryota")))/sum(otu_table(ps))
sum(otu_table(subset_taxa(ps, family%in%"f__Eimeriorina")))/sum(otu_table(ps))

                                        #filter by family

ps_eimf <- subset_taxa(ps, family%in%"f__Eimeriorina")


#ps_eimg <- subset_taxa(ps, genus%in%"g__Eimeria")


sdt18S <- data.frame(sample_sums(otu_table(ps_eimf)))

sdt18S[,2] <- rownames(sdt18S)

sdt18S[,3] <- as.data.frame(sample_sums(ps))

colnames(sdt18S) <- c("Reads_Eimeriori18S", "labels", "Reads_total18S")

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

#source qd-Eim.R to get sdtmulti table.... it's gonna break here because I am lazy

### Join all the data in the same dataframe 
sdt<- join(sample.data, sdt, by="labels") ## First sample data and multimarker read data 
sdt<- join(sdt, data.inf.exp, by="labels") ## then qPCR data
sdt <- join(sdt, sdt18S, by="labels")
sdt <- join(sdt, sdtmulti, by="labels")

sdt <- unique(sdt)

nrow(sdt)

sdt$dpi<- as.factor(sdt$dpi)

opg_multi <-  ggplot(sdt, aes(log(1+OPG), log(1+Reads_Eimeriorina)))+
  geom_jitter(shape=21, position=position_jitter(0.2), size=4, aes(fill= dpi), color= "black", alpha=0.7)+
  labs(tag= "a)")+
  theme_bw()+
  theme(text = element_text(size=16))

opg_qpcr <- ggplot(sdt, aes(log(1+OPG), log(1+Genome_copies_mean)))+
  geom_jitter(shape=21, position=position_jitter(0.2), size=4, aes(fill= dpi), color= "black", alpha=0.7)+
  labs(tag= "A)")+
  theme_bw()+
  theme(text = element_text(size=16))

qpcr_multi <-  ggplot(sdt, aes(log(1+Genome_copies_mean), log(1+Reads_Eimeriorina)))+
    geom_jitter(shape=21, position=position_jitter(0.2), size=4, aes(fill= dpi), color= "black", alpha=0.7)+
        xlab("Genome copies (log 1+)")+
    ylab("multiamplicon Eimeriorina (log1+)")+
    labs(tag= "a)")+
  annotate(geom="text", x=5, y=8.2, label="Spearman rho=0.84, p<0.001")+
  theme_bw()+
  theme(text = element_text(size=16))

cor.test(sdt$Genome_copies_mean, sdt$Reads_Eimeriorina, method="spearman")

cor.test(sdt$Genome_copies_mean, sdt$Reads_Eimeriori18S, method="spearman")

cor.test(sdt$Reads_Eimeriori18S, sdt$Reads_Eimeriorina,method="spearman")



qpcr_18S <-  ggplot(sdt, aes(log(1+Genome_copies_mean), log(1+Reads_Eimeriori18S)))+
    geom_jitter(shape=21, position=position_jitter(0.2), size=4, aes(fill= dpi), color= "black", alpha=0.7)+
    xlab("Genome copies (log 1+)")+
    ylab("single 18S Eimeriorina (log 1+)")+
    labs(tag= "b)")+
    annotate(geom="text", x=5, y=8.2, label="Spearman rho=0.90, p<0.001")+
  theme_bw()+
  theme(text = element_text(size=16))

multi18 <- ggplot(sdt, aes(log(1+Reads_Eimeriori18S), log(1+Reads_Eimeriorina))) +
    geom_jitter(shape=21, position=position_jitter(0.2), size=4, aes(fill=dpi), color="black", alpha=0.7) +
    ylab("multiamplicon Eimeriorina (log 1+)")+
    xlab("single 18S Eimeriorina (log 1+)")+
      annotate(geom="text", x=3, y=8.2, label="Spearman rho=0.86, p<0.001")+
    labs(tag="c") +
    theme_bw() +
    theme(text = element_text(size=16))

multi18

qpcr_18S

opg_multi

qpcr_multi

library(gridExtra)

pdf("18S_multidna.pdf", width=5, height=5)
multi18
dev.off()

pdf("dna_comp.pdf", width=15, height=5)
grid.arrange(qpcr_multi, qpcr_18S, multi18, ncol=3)
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

setwd("/SAN/Susanas_den/gitProj/LabQiime2/single_18S/")

library(phyloseq)
library(ggplot2)
library(ape)


ps_ma<-readRDS(file="/SAN/Susanas_den/gitProj/LabQiime2/single_18S/tmp/PhyloSeqData_ma18S_FullRun.Rds")

ps_18S<-readRDS(file="/SAN/Susanas_den/gitProj/LabQiime2/single_18S/tmp/PhyloSeqData18S.Rds")



ps_ma_eimf <- subset_taxa(ps_ma, family%in%"Eimeriidae")
ps_18S_eimf <- subset_taxa(ps_18S, family%in%"Eimeriidae")
ps_ma_eimg <- subset_taxa(ps_ma, genus%in%"Eimeria")
ps_18S_eimg <- subset_taxa(ps_18S, genus%in%"Eimeria")


sdtma <- data.frame(sample_sums(otu_table(ps_ma_eimf)))
sdtma[,2] <- rownames(sdtma)
sdtma[,3] <- as.data.frame(sample_sums(ps_ma))
colnames(sdtma) <- c("Reads_Eimeriidaema", "labels", "Reads_totalma")

sdt18S <- data.frame(sample_sums(otu_table(ps_18S_eimf)))
sdt18S[,2] <- rownames(sdt18S)
sdt18S[,3] <- as.data.frame(sample_sums(ps_18S))
colnames(sdt18S) <- c("Reads_Eimeriidae18S", "labels", "Reads_total18S")


if(!exists("data.inf.exp")){
  data.inf.exp<- read.csv(file="/SAN/Victors_playground/Eimeria_microbiome/qPCR/sample_data_qPCR.csv")
}

##Keep useful information
##Get unique labels from qPCR data
keeps <- c("labels", "Qty_mean", "Genome_copies_mean","Tm_mean", "Infection")
data.inf.exp <- data.inf.exp[,colnames(data.inf.exp)%in%keeps]

#source qd-Eim.R to get sdtmulti table.... it's gonna break here because I am lazy
### Join all the data in the same dataframe 
sdt <- join(sample.data, sdtma, by="labels") # and single amplicon data
sdt<- join(sdt, data.inf.exp, by="labels") ## join with qPCR data
sdt <- join(sdt, sdt18S, by="labels") # and single amplicon data
sdt <- unique(sdt)

sdt$dpi<- as.factor(sdt$dpi)
cor.test(sdt$Genome_copies_mean, sdt$Reads_Eimeriidaema, method="spearman")
cor.test(sdt$Genome_copies_mean, sdt$Reads_Eimeriidae18S, method="spearman")
cor.test(sdt$Reads_Eimeriidae18S, sdt$Reads_Eimeriidaema,method="spearman")

qpcr_ma <-ggplot(sdt, aes(log(1+Genome_copies_mean), log(1+Reads_Eimeriidaema)))+
    geom_jitter(shape=21, position=position_jitter(0.2), size=4, aes(fill= dpi), color= "black", alpha=0.7)+
    xlab("Genome copies (log 1+)")+
    ylab("multiamplicon Eimeriidae (log1+)")+
    ggtitle("MA processed") +
labs(tag= "a)")+
    annotate(geom="text", x=5, y=2.5, label="Spearman rho=0.12, p=0.06")+
    theme_bw()+
    theme(text = element_text(size=16))

qpcr_18S <-  ggplot(sdt, aes(log(1+Genome_copies_mean), log(1+Reads_Eimeriidae18S)))+
    geom_jitter(shape=21, position=position_jitter(0.2), size=4, aes(fill= dpi), color= "black", alpha=0.7)+
    xlab("Genome copies (log 1+)")+
    ylab("single 18S Eimeriidae (log 1+)")+
    ggtitle("MA processed") +
    labs(tag= "b)")+
    annotate(geom="text", x=5, y=8.2, label="Spearman rho=0.90, p<0.001")+
  theme_bw()+
  theme(text = element_text(size=16))

multi18 <- ggplot(sdt, aes(log(1+Reads_Eimeriidae18S), log(1+Reads_Eimeriidaema))) +
    geom_jitter(shape=21, position=position_jitter(0.2), size=4, aes(fill=dpi), color="black", alpha=0.7) +
    ylab("multiamplicon Eimeriidae (log 1+)")+
    xlab("single 18S Eimeriidae (log 1+)")+
    annotate(geom="text", x=3, y=3, label="Spearman rho=0.15, p=0.03")+
    labs(tag="c)") +
    ggtitle("MA processed") +
    theme_bw() +
    theme(text = element_text(size=16))

library(gridExtra)

pdf("tmp/MA_dna_comp.pdf", width=15, height=5)
grid.arrange(qpcr_ma, qpcr_18S, multi18, ncol=3)
dev.off()

### now with the deblur data

ps_d_18S <- readRDS(file="/SAN/Susanas_den/gitProj/LabQiime2/single_18S/tmp/deblur_phyloseq_single_18S.RDS")
ps_d_18S_eimf <- subset_taxa(ps_d_18S, family%in%"f__Eimeriorina")

sdt_d18S <- data.frame(sample_sums(otu_table(ps_d_18S_eimf)))
sdt_d18S[,2] <- rownames(sdt_d18S)
sdt_d18S[,3] <- as.data.frame(sample_sums(ps_d_18S))
colnames(sdt_d18S) <- c("Reads_Eimeriona_d_18S", "labels", "Reads_total_d_18S")

ps_d <- readRDS(file="/SAN/Susanas_den/gitProj/LabQiime2/merged/tmp/deblur_phyloseq_merged.RDS")
ps_d_eimf <- subset_taxa(ps_d, family%in%"f__Eimeriorina")
ps_d_eimg <- subset_taxa(ps_d, genus%in%"g__Eimeria")


sdtmulti <- data.frame(sample_sums(otu_table(ps_d_eimf)))
sdtmulti[,2] <- rownames(sdtmulti)
sdtmulti[,3] <- as.data.frame(sample_sums(ps_d_eimg))
sdtmulti[,4] <- as.data.frame(sample_sums(ps_d))
colnames(sdtmulti) <- c("Reads_Eimeriorina_d", "labels", "Reads_Eimeria_d", "Reads_total_d")

## and now MA sorted 18S deblur processed
ps_ma_d_18S <- readRDS(file="/SAN/Susanas_den/gitProj/LabQiime2/sorted/tmp/18SmaSorted_deblur_phyloseq.RDS")
ps_ma_d_18S_eimf <- subset_taxa(ps_ma_d_18S, family%in%"f__Eimeriorina")
sdt_ma_d_18S <- data.frame(sample_sums(otu_table(ps_ma_d_18S_eimf)))
sdt_ma_d_18S[,2] <- rownames(sdt_ma_d_18S)
sdt_ma_d_18S[,3] <- as.data.frame(sample_sums(ps_ma_d_18S))
colnames(sdt_ma_d_18S) <- c("Reads_Eimeriidae_mad18S", "labels", "Reads_total_mad18S")

## finally with MA sorted 18S MA processed
# sort amplicons
PS.l <- readRDS(file="/SAN/Susanas_den/gitProj/LabQiime2/single_18S/tmp/PhyloSeqList_ma18S_FullRun.Rds")
PS.l.eimf <- subset_taxa(ps_d_18S, family%in%"f__Eimeriorina")
test <- unlist(PS.l)
ps18S1 <- PS.l$'27M_F_98_F.Klin0341_CR_18_R'
ps18S2 <- PS.l$'515F_Y_118_F.806R_118_R'
ps18S3 <- PS.l$wang1141_13_F.Nem_0425_6_3_R
psMA18S <- merge_phyloseq(ps18S1, ps18S2, ps18S3)

psMA18S

ps_mama_18S_eimf <- subset_taxa(psMA18S, family%in%"Eimeriidae")
ps_mama_18S_eimf
sdt_mama_18S <- data.frame(sample_sums(otu_table(ps_mama_18S_eimf)))
sdt_mama_18S[,2] <- rownames(sdt_mama_18S)
sdt_mama_18S[,3] <- as.data.frame(sample_sums(psMA18S))
colnames(sdt_mama_18S) <- c("Reads_Eimeriidae_mama18S", "labels", "Reads_total_mama18S")


### Join all the data in the same dataframe 

sdt <- join(sdt, sdt_d18S, by="labels") # and single amplicon data

sdt<- join(sdt, sdtmulti, by="labels") ## join with qPCR data

sdt <- join(sdt, sdt_ma_d_18S, by="labels") ## join with ma sorted 18S deblur processed data

sdt <- join(sdt, sdt_mama_18S, by="labels") ## join with ma sorted ma processed 18S data

sdt <- unique(sdt)
head(sdt)

sdt$dpi<- as.factor(sdt$dpi)
cor.test(sdt$Genome_copies_mean, sdt$Reads_Eimeriorina_d, method="spearman")
cor.test(sdt$Genome_copies_mean, sdt$Reads_Eimeriona_d_18S, method="spearman")
cor.test(sdt$Reads_Eimeriona_d_18S, sdt$Reads_Eimeriorina_d,method="spearman")

qpcr_ma_d <-ggplot(sdt, aes(log(1+Genome_copies_mean), log(1+Reads_Eimeriorina_d)))+
    geom_jitter(shape=21, position=position_jitter(0.2), size=4, aes(fill= dpi), color= "black", alpha=0.7)+
    xlab("Genome copies (log 1+)")+
    ylab("multiamplicon Eimeriidae (log1+)")+
    ggtitle("deblur processed") +
labs(tag= "d)")+
    annotate(geom="text", x=5, y=8.2, label="Spearman rho=0.0.84, p<0.001")+
    theme_bw()+
    theme(text = element_text(size=16))

qpcr_18S_d <-  ggplot(sdt, aes(log(1+Genome_copies_mean), log(1+Reads_Eimeriona_d_18S)))+
    geom_jitter(shape=21, position=position_jitter(0.2), size=4, aes(fill= dpi), color= "black", alpha=0.7)+
    xlab("Genome copies (log 1+)")+
    ylab("single 18S Eimeriidae (log 1+)")+
    ggtitle("deblur processed") +
    labs(tag= "e)")+
    annotate(geom="text", x=5, y=8.2, label="Spearman rho=0.90, p<0.001")+
  theme_bw()+
  theme(text = element_text(size=16))

multid <- ggplot(sdt, aes(log(1+Reads_Eimeriona_d_18S), log(1+Reads_Eimeriorina_d))) +
    geom_jitter(shape=21, position=position_jitter(0.2), size=4, aes(fill=dpi), color="black", alpha=0.7) +
    ylab("multiamplicon Eimeriidae (log 1+)")+
    xlab("single 18S Eimeriidae (log 1+)")+
    annotate(geom="text", x=5, y=8.2, label="Spearman rho=0.86, p<0.001")+
    labs(tag="f)") +
    ggtitle("deblur processed") +
    theme_bw() +
    theme(text = element_text(size=16))

head(sdt)

cor.test(sdt$Reads_Eimeriidae18S, sdt$Reads_Eimeriona_d_18S, method="spearman")
cor.test(sdt$Reads_Eimeriidaema, sdt$Reads_Eimeriorina_d,method="spearman")

ma_d_ma <- ggplot(sdt, aes(log(1+Reads_Eimeriidaema), log(1+Reads_Eimeriorina_d))) +
    geom_jitter(shape=21, position=position_jitter(0.2), size=4, aes(fill=dpi), color="black", alpha=0.7) +
    ylab("deblur multiamplicon Eimeriidae (log 1+)")+
    xlab("MA multiamplicon Eimeriidae (log 1+)")+
    annotate(geom="text", x=5, y=8.2, label="Spearman rho=0.21, p=0.001")+
    labs(tag="j)") +
    ggtitle("MA and deblur processed") +
    theme_bw() +
    theme(text = element_text(size=16))

ma_18S_d_ma <- ggplot(sdt, aes(log(1+Reads_Eimeriidae18S), log(1+Reads_Eimeriona_d_18S))) +
    geom_jitter(shape=21, position=position_jitter(0.2), size=4, aes(fill=dpi), color="black", alpha=0.7) +
    ylab("deblur single 18S Eimeriidae (log 1+)")+
    xlab("MA single 18S Eimeriidae (log 1+)")+
    annotate(geom="text", x=5, y=8.2, label="Spearman rho=0.97, p<0.001")+
    labs(tag="k)") +
    ggtitle("MA and deblur processed") +
    theme_bw() +
    theme(text = element_text(size=16))

cor.test(sdt$Reads_Eimeriidae_mama18S, sdt$Reads_Eimeriidae_mad18S, method="spearman")
cor.test(sdt$Genome_copies_mean, sdt$Reads_Eimeriidae_mad18S, method="spearman")
cor.test(sdt$Genome_copies_mean, sdt$Reads_Eimeriidae_mama18S, method="spearman")

qpcr_mad18S <- ggplot(sdt, aes(log(1+Genome_copies_mean), log(1+Reads_Eimeriidae_mad18S))) +
    geom_jitter(shape=21, position=position_jitter(0.2), size=4, aes(fill=dpi), color="black", alpha=0.7) +
    ylab("18S MA + deblur (log 1+)")+
    xlab("qpcr Genome copies (log 1+)")+
    annotate(geom="text", x=5, y=8.2, label="Spearman rho=0.82, p<0.001")+
    labs(tag="g)") +
    ggtitle("multiamplicon MA sorted and deblur processed") +
    theme_bw() +
    theme(text = element_text(size=16))

qpcr_mama18S <- ggplot(sdt, aes(log(1+Genome_copies_mean), log(1+Reads_Eimeriidae_mama18S))) +
    geom_jitter(shape=21, position=position_jitter(0.2), size=4, aes(fill=dpi), color="black", alpha=0.7) +
    ylab("18S sorted and processed with ma (log 1+)")+
    xlab("qpcr Genome copies (log 1+)")+
    annotate(geom="text", x=5, y=8.2, label="Spearman rho=0.82, p<0.001")+
    labs(tag="h)") +
    ggtitle("multiamplicon MA sorted and processed") +
    theme_bw() +
    theme(text = element_text(size=16))

qpcr_mad18S

mama18S_mad18S <- ggplot(sdt, aes(log(1+Reads_Eimeriidae_mama18S), log(1+Reads_Eimeriidae_mad18S))) +
    geom_jitter(shape=21, position=position_jitter(0.2), size=4, aes(fill=dpi), color="black", alpha=0.7) +
    ylab("18S sorted with ma and deblur processed (log 1+)")+
    xlab("18S sorted and processed with ma (log 1+)")+
    annotate(geom="text", x=5, y=8.2, label="Spearman rho=0.82, p<0.001")+
    labs(tag="i)") +
    ggtitle("multiamplicon MA sorted and MA or deblur processed") +
    theme_bw() +
    theme(text = element_text(size=16))

pdf("tmp/MA_d_MAD_dna_comp.pdf", width=20, height=25)
grid.arrange(qpcr_ma, qpcr_18S, multi18, qpcr_ma_d, qpcr_18S_d, multid,
             qpcr_mad18S, qpcr_mama18S, mama18S_mad18S, ma_d_ma, ma_18S_d_ma, ncol=3, nrow=4)
dev.off()

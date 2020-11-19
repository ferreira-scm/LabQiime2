## Please uncomment the first time you run this and re-install packages
##Sequence cleaning and Multiamplicon pipeline for Access Array of Eimeria infection experiment

library(devtools)

devtools::install_github("derele/MultiAmplicon", force= T)

                                        #devtools::install_github("derele/dada2", force= T)

library("lifecycle", lib.loc="/usr/local/lib/R/site-library") 
library("ggplot2")

library("MultiAmplicon")

library("reshape")
library("phyloseq")
library("data.table")
library("taxonomizr")
library("taxize")
library("parallel")

path <- "/SAN/Victors_playground/Eimeria_microbiome/Multimarker/2018_22_Eie_FullRun_1/"

fastqFiles <- list.files(path, pattern=".fastq.gz$", full.names=TRUE) #take all fastaq files from the folder 
fastqF <- grep("_R1_001.fastq.gz", fastqFiles, value = TRUE) #separate the forward reads
fastqR <- grep("_R2_001.fastq.gz", fastqFiles, value = TRUE) #separate the reverse reads 

samples <- gsub("_S\\d+_L001_R1_001.fastq\\.gz", "\\1", basename(fastqF))
samples<- gsub("S\\d+-", "\\1", basename(samples))
samples<-gsub("-", "_", basename(samples))

names(fastqF) <- names(fastqR) <- samples
files <- PairedReadFileSet(fastqF, fastqR)

#Preparation of primer file ### Here stats the Multiamplicon pipeline from Emanuel

#Primers used in the arrays 
ptable <- read.csv(file = "/SAN/Victors_playground/Eimeria_microbiome/Multimarker/primer.file.multi.csv", sep=",", header=TRUE, stringsAsFactors=FALSE)
primerF <- ptable[, "TS.SequenceF"]
primerR <- ptable[, "TS.SequenceR"]
names(primerF) <- as.character(ptable[, "corrected.NameF"])
names(primerR) <- as.character(ptable[, "corrected.NameR"])

primer <- PrimerPairsSet(primerF, primerR)

##Multi amplicon pipeline
MA <- MultiAmplicon(primer, files)

filedir <- "/SAN/Susanas_den/gitProj/LabQiime2/sorted/tmp/stratified_ma_Fullrun"

if(dir.exists(filedir)) unlink(filedir, recursive=TRUE)

MA <- sortAmplicons(MA, n=1e+05, filedir=filedir) ## This step sort the reads into amplicons based on the number of primer pairs

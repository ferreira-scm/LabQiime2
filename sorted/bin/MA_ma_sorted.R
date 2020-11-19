## Please uncomment the first time you run this and re-install packages
##Sequence cleaning and Multiamplicon pipeline for Access Array of Eimeria infection experiment
## require(devtools)
## devtools::install_github("derele/MultiAmplicon", force= T)
## devtools::install_github("derele/dada2", force= T)
library("lifecycle", lib.loc="/usr/local/lib/R/site-library") 
library("ggplot2")
library("MultiAmplicon")
library("reshape")
library("phyloseq")
library("data.table")
library("taxonomizr")
library("taxize")
library("parallel")

## re-run or use pre-computed results for different parts of the pipeline:
## Set to FALSE to use pre-computed and saved results, TRUE to redo analyses.
doFilter <- TRUE
doMultiAmp <- TRUE
doTax <- TRUE
doPhyloseq<- TRUE
## But remember: if you change the MultiAmplicon Analysis, the
## taxonomic annotation might be out of sync...

###################Full run Microbiome#######################
#Preparation of files

##These are the same steps that are followed by the DADA2 pipeline

path <- "/SAN/Victors_playground/Eimeria_microbiome/Multimarker/2018_22_Eie_FullRun_1/"

fastqFiles <- list.files(path, pattern=".fastq.gz$", full.names=TRUE) #take all fastaq files from the folder 
fastqF <- grep("_R1_001.fastq.gz", fastqFiles, value = TRUE) #separate the forward reads
fastqR <- grep("_R2_001.fastq.gz", fastqFiles, value = TRUE) #separate the reverse reads 


samples <- gsub("_S\\d+_L001_R1_001.fastq\\.gz", "\\1", basename(fastqF))
samples<- gsub("S\\d+-", "\\1", basename(samples))
samples<-gsub("-", "_", basename(samples))

#Extra step in the pipeline: quality plots of the reads 
## plotQualityProfile(fastqF[[1]])
## plotQualityProfile(fastqF[[2]])
## plotQualityProfile(fastqR[[1]])
## plotQualityProfile(fastqR[[2]])

#Creation of a folder for filtrated reads 
filt_path <- "/SAN/Susanas_den/gitProj/LabQiime2/single_18S/tmp/filtered_ma18S_FullRun/"
#Pipeline filtration 
if(!file_test("-d", filt_path)) dir.create(filt_path)

filtFs <- file.path(filt_path, paste0(samples, "_F_filt.fastq.gz"))

names(filtFs) <- samples
filtRs <- file.path(filt_path, paste0(samples, "_R_filt.fastq.gz"))
names(filtRs) <- samples

## some files will be filtered out completely, therefore allowing 50
## files less present and still don't redo filtering # full run 150 170 test run 150 150
if(doFilter){
lapply(seq_along(fastqF),  function (i) {
    filterAndTrim(fastqF[i], filtFs[i], fastqR[i], filtRs[i],
                  truncLen=c(200,200), 
                  maxN=0, maxEE=2, truncQ=2, rm.phix = TRUE,
                  compress=TRUE, verbose=TRUE)})
}

names(filtFs) <- names(filtRs) <- samples
files <- PairedReadFileSet(filtFs, filtRs)

#Preparation of primer file ### Here stats the Multiamplicon pipeline from Emanuel

#Primers used in the arrays 
ptable <- read.csv(file = "/SAN/Victors_playground/Eimeria_microbiome/Multimarker/primer.file.multi.csv", sep=",", header=TRUE, stringsAsFactors=FALSE)
primerF <- ptable[, "TS.SequenceF"]
primerR <- ptable[, "TS.SequenceR"]
names(primerF) <- as.character(ptable[, "corrected.NameF"])
names(primerR) <- as.character(ptable[, "corrected.NameR"])

primer <- PrimerPairsSet(primerF, primerR)

##Multi amplicon pipeline
if(doMultiAmp){
  MA <- MultiAmplicon(primer, files)
  filedir <- "/SAN/Susanas_den/gitProj/LabQiime2/single_18S/tmp/stratified_ma18S_Fullrun"
  if(dir.exists(filedir)) unlink(filedir, recursive=TRUE)
  MA <- sortAmplicons(MA, n=1e+05, filedir=filedir) ## This step sort the reads into amplicons based on the number of primer pairs
  errF <-  learnErrors(unlist(getStratifiedFilesF(MA)), nbase=1e8,
                       verbose=0, multithread = 12)
  errR <- learnErrors(unlist(getStratifiedFilesR(MA)), nbase=1e8,
                      verbose=0, multithread = 12)
  MA <- derepMulti(MA, mc.cores=12) 
  MA <- dadaMulti(MA, Ferr=errF, Rerr=errR,  pool=FALSE,
                  verbose=0, mc.cores=12)
  MA <- mergeMulti(MA, mc.cores=12) 
  propMerged <- MultiAmplicon::calcPropMerged(MA)
  MA <- mergeMulti(MA, mc.cores=12) 
  MA <- makeSequenceTableMulti(MA, mc.cores=12) 
  MA <- removeChimeraMulti(MA, mc.cores=12)
  #saveRDS(MA, "/SAN/Victors_playground/Eimeria_microbiome/Multimarker/MA_Multi_TestRun.RDS")
  saveRDS(MA, "/SAN/Susanas_den/gitProj/LabQiime2/single_18S/tmp/MA_ma18S_Fullrun.RDS")
} else{
  #MA <- readRDS("/SAN/Victors_playground/Eimeria_microbiome/MA_Multi_TestRun.RDS")
  MA <- readRDS("/SAN/Susanas_den/gitProj/LabQiime2/single_18S/tmp/MA_ma18S_Fullrun.RDS")
}

getPipelineSummary(MA) 

plotAmpliconNumbers(MA)

###New taxonomic assignment 
if(doTax){
MA <- blastTaxAnnot(MA,
                    db = "/SAN/db/blastdb/nt/nt",
                    negative_gilist = "/SAN/db/blastdb/uncultured.gi",
                    infasta = "/SAN/Susanas_den/gitProj/LabQiime2/single_18S/tmp/in_ma18S_fullrun",
                    outblast = "/SAN/Susanas_den/gitProj/LabQiime2/single_18S/tmp/out_ma18S_fullrun.fasta",
                    taxonSQL = "/SAN/db/taxonomy/taxonomizr.sql", 
                    num_threads = 20) ##Change for use more power!!
saveRDS(MA, file="/SAN/Susanas_den/gitProj/LabQiime2/single_18S/tmp/MA_ma18STax_Fullrun.Rds")
} else {
    MA <- readRDS("/SAN/Susanas_den/gitProj/LabQiime2/single_18S/tmp/MA_ma18STax_Fullrun.Rds")  
}

### Add sample information
if(!exists("sample.data")){
    source("/SAN/Susanas_den/EimeriaMicrobiome/R/1_Data_preparation.R")
}

MA <- addSampleData(MA, sample.data)
#saveRDS(MA2, file="/SAN/Victors_playground/Eimeria_microbiome/Multimarker/MASample_TestRun.Rds")  ###START from here now! 
saveRDS(MA, file="/SAN/Susanas_den/gitProj/LabQiime2/single_18S/tmp/MA_ma18S_SAmple_Fullrun.Rds")
##To phyloseq

if(doPhyloseq){
##Sample data
PS <- toPhyloseq(MA, colnames(MA))
sum(otu_table(PS)) ##Total denoised reads = 3,610,637 for FullRun
##Primer data
## just sorting out primers whithout any taxannot
MA <- MA[which( !unlist(lapply(MA@taxonTable, is.null))), ] ##Make the next function work 
PS.l <- toPhyloseq(MA, colnames(MA),  multi2Single=FALSE) 
saveRDS(PS.l, file="/SAN/Susanas_den/gitProj/LabQiime2/single_18S/tmp/PhyloSeqList_ma18S_FullRun.Rds")
saveRDS(PS, file="/SAN/Susanas_den/gitProj/LabQiime2/single_18S/tmp/PhyloSeqData_ma18S_FullRun.Rds") ###For Sample analysis (Susana and Victor)
}else{
####Merge phyloseq objects   
    PS <-  readRDS(file="/SAN/Susanas_den/gitProj/LabQiime2/single_18S/tmp/PhyloSeqData_ma18S_FullRun.Rds")
}

# sort amplicons
PS.l <- readRDS(file="/SAN/Susanas_den/gitProj/LabQiime2/single_18S/tmp/PhyloSeqList_ma18S_FullRun.Rds")

PS.l.eimf <- subset_taxa(ps_d_18S, family%in%"f__Eimeriorina")

test <- unlist(PS.l)

ps18S1 <- PS.l$'27M_F_98_F.Klin0341_CR_18_R'

ps18S2 <- PS.l$'515F_Y_118_F.806R_118_R'

ps18S3 <- PS.l$wang1141_13_F.Nem_0425_6_3_R

psMA18S <- merge_phyloseq(ps18S1, ps18S2, ps18S3)

psMA18S

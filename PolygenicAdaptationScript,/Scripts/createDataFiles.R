###############################################################################
### Use this script to run the Berg and Coop method on both_conditions      ###
### height data.                                                            ###
###############################################################################

### Timothy M. Beissinger
### 10-21-2014

###############################################################################
### Create freqs file                                                       ###
###############################################################################
### Run loop to load highlow SNPs by chromosome
fulldata <- read.table("../Data/highlow_gbs_hapmaps 2/C08JYACXX_2_c1.hmp.txt",header=T,comment.char="",stringsAsFactors=F) #read chromosome 1
fulldata$chr <- 1
for(chr in 2:10){
    print(chr)
    chrXhighlow <- read.table(paste("../Data/highlow_gbs_hapmaps 2/C08JYACXX_2_c",chr,".hmp.txt",sep=""),header=T,comment.char="",stringsAsFactors=F)
    chrXhighlow$chr <- chr
    fulldata <- rbind(fulldata,chrXhighlow)
 }

# make backup data frame
fulldataFrame <- fulldata

### Replace hapmap code with nucleotides...
fulldata <- as.matrix(fulldataFrame)
fulldata[which(fulldata=="N")] <- NA
fulldata[which(fulldata=="X")] <- NA
fulldata[which(fulldata=="A")] <- "A/A"
fulldata[which(fulldata=="C")] <- "C/C"
fulldata[which(fulldata=="G")] <- "G/G"
fulldata[which(fulldata=="T")] <- "T/T"
fulldata[which(fulldata=="K")] <- "G/T"
fulldata[which(fulldata=="M")] <- "A/C"
fulldata[which(fulldata=="R")] <- "A/G"
fulldata[which(fulldata=="S")] <- "C/G"
fulldata[which(fulldata=="W")] <- "A/T"
fulldata[which(fulldata=="Y")] <- "C/T"

### make fulldata a data frame again
fulldata <- as.data.frame(fulldata,stringsAsFactors=F)

### Load membership information, created BY HAND!!
load("../Data/highlowMembership.Robj")

### Make membership vector
membership <- c("NA","NA","NA","NA",highlowMembership[,2],"NA") # perfectly corresponds to fulldata

### Determine allele frequencies
freqs <- matrix(NA,nrow=nrow(fulldata),ncol=12)
colnames(freqs) <- c("rs","Alleles","Chr","pos","MHcount","MHfreq","MLcount","MLfreq","SAHcount","SAHfreq","SALcount","SALfreq")
freqs <- as.data.frame(freqs)

for(i in 1:nrow(fulldata)){
    print(i)
    freqs$rs[i] <- fulldata$rs[i]
    freqs$Alleles[i] <- fulldata$alleles[i]
    freqs$Chr[i] <- as.numeric(fulldata$chr[i])
    freqs$pos[i] <- as.numeric(fulldata$pos[i])
    allele <- substr(fulldata$alleles[i],1,1)

    MHalleles <- unlist(strsplit(as.character(fulldata[i,which(membership=="MH")]),split="/"))
    freqs$MHcount[i] <- length(which(is.na(MHalleles)==F))
    freqs$MHfreq[i] <- length(which(MHalleles==allele))/length(which(is.na(MHalleles)==F))

    MLalleles <- unlist(strsplit(as.character(fulldata[i,which(membership=="ML")]),split="/"))
    freqs$MLcount[i] <- length(which(is.na(MLalleles)==F))
    freqs$MLfreq[i] <- length(which(MLalleles==allele))/length(which(is.na(MLalleles)==F))

    SAHalleles <- unlist(strsplit(as.character(fulldata[i,which(membership=="SAH")]),split="/"))
    freqs$SAHcount[i] <- length(which(is.na(SAHalleles)==F))
    freqs$SAHfreq[i] <- length(which(SAHalleles==allele))/length(which(is.na(SAHalleles)==F))

    SALalleles <- unlist(strsplit(as.character(fulldata[i,which(membership=="SAL")]),split="/"))
    freqs$SALcount[i] <- length(which(is.na(SALalleles)==F))
    freqs$SALfreq[i] <- length(which(SALalleles==allele))/length(which(is.na(SALalleles)==F))
}

### Backup freqs file
freqsFull <- freqs


### Find columns in freqs file with NAs
NAcols <- which(is.na(freqs$MHfreq) | is.na(freqs$MLfreq) | is.na(freqs$SALfreq) | is.na(freqs$SAHfreq))
freqs <- freqs[-NAcols,]

### Write freqs file as Robject for analysis with other software
save(freqs,file="highlow-freqs.Robj")

### CHECKPOINT ###
save.image("createDataFiles.RData")

#### Create freqs file
freqs.file <- matrix(NA,nrow=4*nrow(freqs),ncol=8)
freqs.file <- as.data.frame(freqs.file)
names(freqs.file) <- c("SNP","CLST","A1","A2","FRQ","POS","CHR","COUNT")

chr <- unlist(strsplit(freqs$rs,split="_"))[seq(1,2*nrow(freqs),2)] #chromosome
chr <- as.numeric(substr(chr,2,10)) #chromosome

pos <- unlist(strsplit(freqs$rs,split="_"))[seq(2,2*nrow(freqs),2)] #position
pos <- as.numeric(pos) #position

freqs.file[1:nrow(freqs),] <- cbind(freqs$rs,rep("MH",nrow(freqs)),substr(freqs$Alleles,1,1),substr(freqs$Alleles,3,3),freqs$MHfreq,pos,chr,freqs$MHcount)
freqs.file[{nrow(freqs)+1}:{2*nrow(freqs)},] <- cbind(freqs$rs,rep("ML",nrow(freqs)),substr(freqs$Alleles,1,1),substr(freqs$Alleles,3,3),freqs$MLfreq,pos,chr,freqs$MLcount)
freqs.file[{2*nrow(freqs)+1}:{3*nrow(freqs)},] <- cbind(freqs$rs,rep("SAH",nrow(freqs)),substr(freqs$Alleles,1,1),substr(freqs$Alleles,3,3),freqs$SAHfreq,pos,chr,freqs$SAHcount)
freqs.file[{3*nrow(freqs)+1}:{4*nrow(freqs)},] <- cbind(freqs$rs,rep("SAL",nrow(freqs)),substr(freqs$Alleles,1,1),substr(freqs$Alleles,3,3),freqs$SALfreq,pos,chr,freqs$SALcount)

freqs.file$CHR <- as.numeric(freqs.file$CHR) #make numeric
freqs.file$POS <- as.numeric(freqs.file$POS) #make numeric
freqs.file$FRQ <- as.numeric(freqs.file$FRQ) #make numeric
freqs.file$COUNT <- as.numeric(freqs.file$COUNT) #make numeric

freqs.file <- freqs.file[order(freqs.file$CHR,freqs.file$POS,freqs.file$CLST),] #put in order, too many SNPs, trim later


###############################################################################
### First, define a function to create gwas data file and freqs file        ###
###############################################################################
makeDataFiles <- function(phenoPath,freqsFile = freqs.file,label){
    print(phenoPath)
    gwasdata <- read.table(phenoPath,sep="\t",header=T,stringsAsFactors=F) #read
    gwasdataSig <- gwasdata[which(gwasdata$pvalCP3 <= 1e-4),] #say significance threshold is 1e-4
    gwasdataSig$FRQ <- gwasdataSig$CountAll2 / {gwasdataSig$CountRef + gwasdataSig$CountAll2}
    gwasDataFile <- gwasdataSig[,c(1,9,11)]
    names(gwasDataFile) <- c("SNP","EFF","FRQ")

    gwasDataFileTrimmed <- gwasDataFile[which(gwasDataFile$SNP %in% freqsFile$SNP),]
    freqsFileTrimmed <- freqsFile[which(freqsFile$SNP %in% gwasDataFileTrimmed$SNP),]
    print(paste("GWAS HITS = ", nrow(gwasDataFileTrimmed),sep=""))

    print("Writing Trimmed Files")
    write.table(gwasDataFileTrimmed,file=paste("Trait_Data/gwas.data.",label,".txt",sep=""),row.names=F,col.names=T,quote=F) # write gwas.data.file
    write.table(freqsFileTrimmed,file=paste("Trait_Data/freqs.file.",label,".txt",sep=""),row.names=F,col.names=T,quote=F) # write freqs.file
}


makeDataFiles100 <- function(phenoPath,freqsFile = freqs.file,label){
    print(phenoPath)
    gwasdata <- read.table(phenoPath,sep="\t",header=T,stringsAsFactors=F) #read
    gwasdataSig <- gwasdata[order(gwasdata$pvalCP3)[1:100],] #Take top 100 as significant
    gwasdataSig$FRQ <- gwasdataSig$CountAll2 / {gwasdataSig$CountRef + gwasdataSig$CountAll2}
    gwasDataFile <- gwasdataSig[,c(1,9,11)]
    names(gwasDataFile) <- c("SNP","EFF","FRQ")

    gwasDataFileTrimmed <- gwasDataFile[which(gwasDataFile$SNP %in% freqsFile$SNP),]
    freqsFileTrimmed <- freqsFile[which(freqsFile$SNP %in% gwasDataFileTrimmed$SNP),]
    print(paste("GWAS HITS = ", nrow(gwasDataFileTrimmed),sep=""))

    print("Writing Trimmed Files")
    write.table(gwasDataFileTrimmed,file=paste("Trait_Data100/gwas.data.",label,".txt",sep=""),row.names=F,col.names=T,quote=F) # write gwas.data.file
    write.table(freqsFileTrimmed,file=paste("Trait_Data100/freqs.file.",label,".txt",sep=""),row.names=F,col.names=T,quote=F) # write freqs.file
}



####################################################
### TRIM AND WRITE freqs.file and gwas.data.file ###
####################################################
for(i in 1:length(dir("GWASResults"))){
makeDataFiles(paste("GWASResults/",dir("GWASResults")[i],sep=""),freqsFile = freqs.file,label=strsplit(dir("GWASResults")[i],split=".txt")[[1]])
}


for(i in 1:length(dir("GWASResults"))){
makeDataFiles100(paste("GWASResults/",dir("GWASResults")[i],sep=""),freqsFile = freqs.file,label=strsplit(dir("GWASResults")[i],split=".txt")[[1]])
}

####################################################
### Write full.dataset.file                      ###
####################################################
write.table(freqs.file,file="Genome_Data/full.dataset.file.txt",row.names=F,col.names=T,quote=F) # write freqs.file


###############################################################################
### Next, create match.pop file                                             ###
###############################################################################

### Run loop to load maize282 SNPs by chromosome
fulldata282 <- read.table("../Data/282_unique_taxa/chr_1.txt",header=T,comment.char="",stringsAsFactors=F) #read chromosome 1
chr1highlow <- freqs[which(freqs$Chr == 1),]#only SNPs with highlow freqs, in "freqs"
fulldata282 <- fulldata282[which(fulldata282$pos %in% chr1highlow$pos),]#only SNPs with highlow freqs

for(chr in 2:10){
    print(chr)
    chrX282 <- read.table(paste("../Data/282_unique_taxa/chr_",chr,".txt",sep=""),header=T,comment.char="",stringsAsFactors=F)
    chrXhighlow <- freqs[which(freqs$Chr == chr),]#only SNPs with highlow freqs, in "freqs"
    chrX282 <- chrX282[which(chrX282$pos %in% chrXhighlow$pos),]#only SNPs with highlow freqs
    fulldata282 <- rbind(fulldata282,chrX282)
 }


# make backup data frame
fulldata282Frame <- fulldata282

### Replace hapmap code with nucleotides...
fulldata282 <- as.matrix(fulldata282Frame)
fulldata282[which(fulldata282=="N")] <- NA
fulldata282[which(fulldata282=="X")] <- NA
fulldata282[which(fulldata282=="A")] <- "A/A"
fulldata282[which(fulldata282=="C")] <- "C/C"
fulldata282[which(fulldata282=="G")] <- "G/G"
fulldata282[which(fulldata282=="T")] <- "T/T"
fulldata282[which(fulldata282=="K")] <- "G/T"
fulldata282[which(fulldata282=="M")] <- "A/C"
fulldata282[which(fulldata282=="R")] <- "A/G"
fulldata282[which(fulldata282=="S")] <- "C/G"
fulldata282[which(fulldata282=="W")] <- "A/T"
fulldata282[which(fulldata282=="Y")] <- "C/T"

### CHECKPOINT ###
save.image("createDataFiles.RData")

### make fulldata282 a data frame again
fulldata282 <- as.data.frame(fulldata282,stringsAsFactors=F)

### Determine allele to track
SNP <- paste("S",as.numeric(fulldata282$chrom),"_",as.numeric(fulldata282$pos),sep="") #name SNPs
freqsSub <- freqs[which(freqs$rs %in% SNP),] # isolate a subset of tracked SNPs
Allele1 <- substr(freqsSub$Alleles,1,1) #record highlow tracked allele

### Determine allele frequencies
freqs282 <- matrix(NA,nrow=nrow(fulldata282),ncol=7)
colnames(freqs282) <- c("Alleles","A1","A2","Chr","pos","count","freq")
freqs282 <- as.data.frame(freqs282)
freqs282$Alleles <- fulldata282$alleles
freqs282$A1 <- Allele1 #from above, corresponds to highlow tracked allele
freqs282$Chr <- fulldata282$chrom
freqs282$pos <- fulldata282$pos

AlleleList <- strsplit(freqs282$Alleles,split="/")


for(i in 1:nrow(fulldata282)){
    print(i)
    allele <- freqs282$A1[i]
    freqs282$A2[i] <- AlleleList[[i]][which(AlleleList[[i]] != allele)]

    locus <- fulldata282[i,11:ncol(fulldata282)]
    alleles <- unlist(strsplit(as.character(locus),split="/"))
    freqs282$count[i] <- length(which(is.na(alleles)==F))
    freqs282$freq[i] <- length(which(alleles==allele))/length(which(is.na(alleles)==F))
}

### CHECKPOINT ###
save.image("createDataFiles.RData")

### Find SNPs for which highlow SNPs are wholly different from 282 SNPs
rem <- c()
for(i in 1:nrow(freqs282)){
    print(i)
    if(length(which(AlleleList[[i]] != freqs282$A1[i])) > 1) rem <- c(rem,i)
}
### Probably no need to remove them.

### Create match.pop.file
match.pop.file <- matrix(NA,nrow=nrow(freqs282),ncol=7)
match.pop.file <- as.data.frame(match.pop.file)
names(match.pop.file) <- c("SNP","CLST","A1","A2","FRQ","POS","CHR")
match.pop.file$SNP <- SNP
match.pop.file$CLST <- "Ames282" #name population
match.pop.file$A1 <- freqs282$A1
match.pop.file$A2 <- freqs282$A2
match.pop.file$FRQ <- freqs282$freq
match.pop.file$POS <- as.numeric(freqs282$pos)
match.pop.file$CHR <- as.numeric(freqs282$Chr)

### Write match.pop.file
write.table(match.pop.file,file="Genome_Data/match.pop.file.txt",row.names=F,col.names=T,quote=F)

save.image("createDataFiles.RData")

rm(list=ls())
source("../Scripts/functions_modified.R")

PolygenicAdaptationFunction ( 
  gwas.data.file = "../../qxfiles/1kwt/gwas_file2.csv", 
  freqs.file = "../../qxfiles/1kwt/freqfile2.csv" , 
  env.var.data.files = list ("../EnvVar/Environment.txt" ) , 
  # Note: you can supply as many env.var.data.files concurrently as you want. If only supplying one file it should still be included in a list, e.g. env.var.data.files = list ( "Example/EnvVar/HGDP_LATS_GLOBAL")
  match.pop.file = "../../qxfiles/1kwt/match_pop_file2.csv",
  full.dataset.file = "../../qxfiles/1kwt/full_dataset_file2.csv" , 
  path = "../../qxfiles/1kwt/output/" , 
  match.categories = c ("FRQ") ,
  match.bins = list ( seq ( 0 , 1, 0.5 )) ,
  cov.SNPs.per.cycle = 5000 , 
  cov.cycles = 1 , 
  null.phenos.per.cycle = 1000 , 
  null.cycles = 1 ,
  load.cov.mat = F ,
  sim.null = T ,
  check.allele.orientation = T
)
######################### DO ANALYSIS FOR GWAS P-VALUES < 1e-4 #########################
Out <- list()
for(i in 1:22){
traitFile = paste("Trait_Data/",dir("Trait_Data")[i+22],sep="")
freqsFile = paste("Trait_Data/",dir("Trait_Data")[i],sep="")
pathDir = strsplit(traitFile,split="[.]")[[1]][3]

Out[[i]] <- PolygenicAdaptationFunction (	gwas.data.file = traitFile,
				freqs.file = freqsFile ,
				env.var.data.files = list ( "EnvVar/Environment.txt") ,
							match.pop.file = "Genome_Data/match.pop.file.txt" ,
							full.dataset.file = "Genome_Data/full.dataset.file.txt" ,
							path = paste("OUTPUT/",pathDir,sep="") ,
							match.categories = c ( "FRQ") ,
							match.bins = list ( seq ( 0 , 1 , 0.02 ) ) ,
							cov.SNPs.per.cycle = 5000 ,
							cov.cycles = 1 ,
							null.phenos.per.cycle = 1000 ,
							null.cycles = 10 ,
							load.cov.mat = T ,
							sim.null = F ,
							check.allele.orientation = F
							)


### Make table of asymptotic Qx p-values
QxPvalues <- matrix(NA,ncol=6,nrow=22)
QxPvalues <- data.frame(QxPvalues)
names(QxPvalues) <- c("Phenotype","Overdispersion-P","MH-P","ML-P","SAH-P","SAL-P")
for(i in 1:22){
QxPvalues[i,1] <- names(Out)[i]
QxPvalues[i,2] <- Out[[i]]$asymptotic.p.vals$Qx
QxPvalues[i,3:6] <- Out[[i]]$asymptotic.p.vals$ind.Z
}


### Plot results
pdf("Summary_Plot_10-4.pdf",height=8.5,width=11)
par(mar=c(7,4,4,2))
plot(-log10(QxPvalues[,2]),xaxt="n",xlab="",ylab=expression("-log"[10]*" p-value"),pch=19,ylim=c(0,max(-log10(QxPvalues[,2:6]))),cex=2,main="Test for polygenic adaptation")
axis(1,at=1:22,labels=F)
text(x=1:22,y=-0.3,labels=QxPvalues[,1],srt=45,adj=1,xpd=T,cex=0.75)
abline(h=-log10(0.05),lty=2)
points(-log10(QxPvalues[,3]),col="red",pch=16)
points(-log10(QxPvalues[,4]),col="orange",pch=16)
points(-log10(QxPvalues[,5]),col="blue",pch=16)
points(-log10(QxPvalues[,6]),col="purple",pch=16)

legend("topright","(x,y)",c("Overdispersion", "MH","ML","SAH","SAL"),pch=19,pt.cex=c(2,1,1,1,1),col=c("black","red","orange","blue","purple"))

dev.off()













######################### DO ANALYSIS FOR LOWEST 100 GWAS P-VALUES #########################
Out100 <- list()
for(i in 1:22){
traitFile = paste("Trait_Data100/",dir("Trait_Data")[i+22],sep="")
freqsFile = paste("Trait_Data100/",dir("Trait_Data")[i],sep="")
pathDir = strsplit(traitFile,split="[.]")[[1]][3]

Out100[[i]] <- PolygenicAdaptationFunction (	gwas.data.file = traitFile,
				freqs.file = freqsFile ,
				env.var.data.files = list ( "EnvVar/Environment.txt") ,
							match.pop.file = "Genome_Data/match.pop.file.txt" ,
							full.dataset.file = "Genome_Data/full.dataset.file.txt" ,
							path = paste("OUTPUT100/",pathDir,sep="") ,
							match.categories = c ( "FRQ") ,
							match.bins = list ( seq ( 0 , 1 , 0.02 ) ) ,
							cov.SNPs.per.cycle = 5000 ,
							cov.cycles = 1 ,
							null.phenos.per.cycle = 1000 ,
							null.cycles = 10 ,
							load.cov.mat = T ,
							sim.null = F ,
							check.allele.orientation = F
							)
names(Out100)[i] <- pathDir
print("------------------------------------------------------------------")
print(Out100[[i]])
}


### Make table of asymptotic Qx p-values
QxPvalues100 <- matrix(NA,ncol=6,nrow=22)
QxPvalues100 <- data.frame(QxPvalues100)
names(QxPvalues100) <- c("Phenotype","Overdispersion-P","MH-P","ML-P","SAH-P","SAL-P")
for(i in 1:22){
QxPvalues100[i,1] <- names(Out100)[i]
QxPvalues100[i,2] <- Out100[[i]]$asymptotic.p.vals$Qx
QxPvalues100[i,3:6] <- Out100[[i]]$asymptotic.p.vals$ind.Z
}


### Plot results
pdf("Summary_Plot_100.pdf",height=8.5,width=11)
par(mar=c(7,4,4,2))
plot(-log10(QxPvalues100[,2]),xaxt="n",xlab="",ylab=expression("-log"[10]*" p-value"),pch=19,ylim=c(0,max(-log10(QxPvalues[,2:6]))),cex=2,main="Test for polygenic adaptation")
axis(1,at=1:22,labels=F)
text(x=1:22,y=-0.3,labels=QxPvalues[,1],srt=45,adj=1,xpd=T,cex=0.75)
abline(h=-log10(0.05),lty=2)
points(-log10(QxPvalues100[,3]),col="red",pch=16)
points(-log10(QxPvalues100[,4]),col="orange",pch=16)
points(-log10(QxPvalues100[,5]),col="blue",pch=16)
points(-log10(QxPvalues100[,6]),col="purple",pch=16)

legend("topright","(x,y)",c("Overdispersion", "MH","ML","SAH","SAL"),pch=19,pt.cex=c(2,1,1,1,1),col=c("black","red","orange","blue","purple"))

dev.off()




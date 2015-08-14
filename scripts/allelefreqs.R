#!/usr/bin/Rscript

rm(list=ls())
load_obj <- function (f) {
  env <- new.env()
  nm <- load(f, env)[1]
  return(env[[nm]])
}

# allelefreqs <- function (filename) {
#   myobj <- load_obj(filename)
#   return(myobj)
# }
# 
# getfreq <- function (geno) {
# 	return(sum(geno) / (length(geno) *2 ))
# }
# 
# calc_freqs <- function (filename) {
# 
# 	alldata <- load_obj(filename)
# 	mycols <- alldata[,-1]
# 	freqs <- apply(mycols, 2, getfreq)
# 	return(freqs)
# }

#running into memory problems running the whole thing so loop by chromosome
allnames = list()
for (x in 1:10) {
	filename <- paste("../data/Namfreqs/cut2uniqueAD_NAM_KIDS_chr", x, ".raw.RData", sep='')
	myobj <- load_obj(filename)
  mynames = names(myobj)
  allnames[[x]] = mynames
}
rm(myobj)
rm(mynames)
fullnames = do.call(c, allnames)
write(fullnames, file = "../data/Namfreqs/GWAS_SNPs.txt")

rm(list=ls())

load_obj <- function (f) {
  env <- new.env()
  nm <- load(f, env)[1]
  return(as.matrix(env[[nm]]))
}

allelefreqs <- function (filename) {
  myobj <- load_obj(filename)
  return(myobj)
}
  
filenames = mapply(paste, "../data/Namfreqs/cut2uniqueAD_NAM_KIDS_chr", 1:10, ".raw.RData", sep='')
allchroms = lapply(files, load_obj)


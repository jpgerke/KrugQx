require(stringr)
require(dplyr)
library(reshape2)
rm(list=ls())

#load the phenotype data
load("../data/1kweight/SNP_effects_ordered.RData")
str(df)

##### load the genome files #####
#for the populations of interest
pop_genome_file <- read.csv("../data/krug_olap_biallelic.csv", as.is= T)
str(pop_genome_file)

#for the gwas panel
gwas_genome_file <- read.csv("../data/nam_olap_krug_biallelelic.csv", as.is = T)
str(gwas_genome_file)

#create an ID for the SNP effects
with_ID <- rowwise(df) %>%
             mutate(SNP = str_c(str_split(term, pattern = "_")[[1]][1:2], collapse = "_"), 
                    Allele = str_split(term, pattern = "_")[[1]][3])

##### match.pop.file
#Format: SNP, CLST, A1, A2, FRQ, POS, CHR, etc.
match_pop_file <- anti_join(gwas_genome_file, with_ID, by="SNP")

#once ready we will want to make note of which allele is the reference allele

#this will be the full.dataset.file
full_dataset_file_shortform <- anti_join(pop_genome_file, with_ID, by="SNP")
longform <- melt(full_dataset_file_shortform, 
                 measure.vars = c("KC0Allele1", "KLSAllele1", "KSSAllele1"), 
                 variable.name = "CLST", 
                 value.name = "FRQ") %>%
            rowwise() %>%
            mutate(CHR = str_split(SNP, pattern = "_")[[1]][1], POS = str_split(SNP, pattern = "_")[[1]][2]) %>%
            rowwise() %>%
            mutate(CHR = str_replace(CHR, pattern='S', replacement = '')[1]) %>%
            select(SNP,CLST,A1,A2,FRQ,POS,CHR)

#####GWAS Hits
#Format:
# SNP A1 A2 EFF FRQ
gwas_data_file <- inner_join(with_ID, gwas_genome_file, by="SNP")
stopifnot( sum(gwas_data_file$Allele == gwas_data_file$A1) == nrow(gwas_data_file))
gwas_data_file_final = select(gwas_data_file, SNP, A1, A2, EFF=estimate, FRQ)

#####FREQS FILE:
#GWAS hits in pops of interest
#Format: SNP, CLST, A1, A2, FRQ, IMP, POS, CHR, etc.
freqs_file <- inner_join(with_ID, pop_genome_file, by="SNP") %>%
  melt(.,
       measure.vars = c("KC0Allele1", "KLSAllele1", "KSSAllele1"), 
       variable.name = "CLST", 
       value.name = "FRQ") %>%
  rowwise() %>%
  mutate(CHR = str_split(SNP, pattern = "_")[[1]][1], POS = str_split(SNP, pattern = "_")[[1]][2]) %>%
  rowwise() %>%
  mutate(CHR = str_replace(CHR, pattern='S', replacement = '')[1]) %>%
  select(SNP, CLST, A1, A2, FRQ, POS, CHR)

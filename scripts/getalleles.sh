#!/bin/bash -l

#SBATCH -J getnam
#SBATCH -e ../logs/getnam-%j-log.err
#SBATCH -o ../logs/getnam-%j-log.out

zcat /group/jrigrp4/Justin_Kate/GBS2.7/sorted_NAM_children.vcf.gz | awk '{print $1, $2, $3}' > ../data/Namfreqs/alleleconversion.txt

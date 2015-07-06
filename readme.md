# KrugQx

### Files needed for Qx

Qx runs on four main genotype input files. 

#### Two of these revolve around the GWAS hits. 

* gwas.data.file: GWAS effect estimates from the GWAS panel. Format:
		
		Name, Allele1, Allele2, Effect, Freq
* freqs.file: allele frequencies and other attributes of the GWAS hits in the test populations. Format:

		Name, CLST, A1, A2, Freq, vars
		
	where `CLST` denotes the population and `vars` are variables to match for the null SNPs. Tim just uses frequency. For the NAM data, we could also use whether the positive effect of major allele came from B73.
	
#### Two are the 'genome' files for all SNPs.

* full.dataset.file: all allele frequencies in all populations, similar in format to freqs.file
* match.pop.file:  the frequencies and matching attributes of all SNPs in the GWAS panel.

#### Qx also uses an environmental variables file.

Of the format `CLUST,ENV,REG`. `REG` corresponds to geographic clusters and can be ignored as we aren't interested in a regional Z score. 

Notes:

* The allele frequencies are the frequencies of allele A1
* The effects are the substitution effects of replacing A2 with A1.

### Available Files 
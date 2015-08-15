
import pickle as pkl

freqs = pkl.load(open("../data/Namfreqs/namfreqdict.pkl", 'rb'))

snpfile = open("../data/Namfreqs/GWAS_SNPS.txt", 'r')

outfile = open('nam_genome_file.csv', 'w')
outfile.write("SNP,A1,A2,FRQ\n")

count = 0
found = 0
for line in snpfile:
	entry = line.rstrip().split('_')
	if entry[0] == "FID":
		continue
	count += 1
	SNP = entry[0] + "_" + entry[1]
	#allele 1 needs to be the allele the effect was calculated for
	myallele = entry[2]
	if SNP in freqs:
		myfreq = freqs[SNP][myallele]
		alleles = freqs[SNP].keys()
		if alleles[0] == myallele:
			A1 = alleles[0]
			A2 = alleles[1]
		else:
			assert alleles[1] == myallele
			A1 = alleles[1]
			A2 = alleles[0]
	#	print SNP, freqs[SNP], myfreq, myallele, A1, A2
		found += 1
		outfile.write(','.join(map(str, [SNP, A1, A2, myfreq])) + '\n')

print count, found

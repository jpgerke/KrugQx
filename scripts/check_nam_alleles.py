
import pickle

alleledict = pickle.load(open("../data/Namfreqs/alleledict.pkl", 'rb'))

outfile = open("../data/Namfreqs/NAM_allfreqs.csv", 'w')
outfile.write("SNP,A1,A2,FRQ\n")

for chrom in xrange(1,11):
	#loop through the frequency file
	infile = open("../data/Namfreqs/frequencies_chr" + str(chrom) + ".csv", 'r')
	print chrom
	for line in infile:
		entry = line.rstrip().split(',')
		#construct the SNP id
		ids = entry[0].split('_')
		snp_id = '_'.join([ids[0], ids[1]])
		#get the two alleles from the dict
		myalleles = alleledict[snp_id]
		#convert dots to dashes
		if ids[2] == '.':
			ids[2] = '-'
		#allele 1 needs to be the allele in ids[2]
		# if neither allele matches the underscore we have a problem
		if ids[2] == myalleles[0]:
			allele1 = myalleles[0]
			allele2 = myalleles[1]

		elif ids[2] == myalleles[1]:
			allele1 = myalleles[1]
			allele2 = myalleles[0]
		else:	
			raise ValueError, "Warning: allele does not match"
		#print out non-indel alleles
		if allele1 != '-' and allele2 != '-':
			outfile.write(','.join([snp_id,allele1,allele2,entry[1]]) + '\n')

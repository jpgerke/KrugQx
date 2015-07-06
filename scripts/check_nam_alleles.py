
import pickle

alleledict = pickle.load(open("../data/Namfreqs/alleledict.pkl", 'rb'))

outfile = open("../data/Namfreqs/allfreqs.csv", 'w')
outfile.write("SNP,A1,A2,FRQ\n")

for chrom in xrange(1,11):
	infile = open("../data/Namfreqs/frequencies_chr" + str(chrom) + ".csv", 'r')

	for line in infile:
		entry = line.rstrip().split(',')
		ids = entry[0].split('_')
		snp_id = '_'.join([ids[0], ids[1]])
		myalleles = alleledict[snp_id]
		if ids[2] == '.':
			ids[2] = '-'
		#assert ids[2] == myalleles[1]
		if ids[2] != myalleles[1]:
				print myalleles, entry
		outfile.write(','.join([snp_id,myalleles[0],myalleles[1],entry[1]]) + '\n')

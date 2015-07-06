
import pickle as pkl

infile = open("../data/Namfreqs/alleles.txt", 'r')

alleledict = {}

#skip the first line
infile.readline()

for line in infile:
	entry = line.rstrip().split()
	#format will be major, minor
	alleledict[entry[0]] = (entry[1], entry[2])

#for key,value in alleledict.iteritems():
#	print key, value
#	break

pkl.dump(alleledict, open("../data/Namfreqs/alleledict.pkl", 'wb'))

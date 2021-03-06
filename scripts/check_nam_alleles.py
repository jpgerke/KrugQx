
import pickle

#construct the dict to convert between v2 and v3
condict = {}
confile = open("../data/Namfreqs/alleleconversion.txt", 'r')
for line in confile:
	if line[0] != '#':
		entry = line.rstrip().split()
		v3 = entry[0] + '_' + entry[1]
		condict[v3] = entry[2]
#print condict

#load the file and skip the header
freqfile = open("../data/Namfreqs/NAMkids_frequencies.frq", 'r')
header = freqfile.readline()

#define the chroms and bases of interest
chroms = map(str, range(1,11))
bases = set(['A', 'C', 'T', 'G'])

#initialize the dictionary to store the data
alleledict = {}

#get down to business
for line in freqfile:
	entry = line.rstrip().split()
	#has to be on a chrom
	if entry[0] not in chroms:
		continue
	#has to be biallelic
	elif int(entry[2]) != 2:
		continue
	#if it is biallelic and on a chrom, proceed:
	else:
		#grab the alleles and their frequencies
		alleles = entry[4::]
		freqs = [x.split(":") for x in alleles]
		
		#sanity check
		assert len(freqs) == 2
		
		#make sure we're dealing with SNPs
		gtype = set([x[0] for x in freqs])
		if gtype.issubset(bases):
			#build a dict for the allele
			v3SNP = entry[0] + "_" + entry[1]
			v2SNP = condict[v3SNP]
			alleledict[v2SNP] = {}
			for item in freqs:
				alleledict[v2SNP][item[0]] = float(item[1])

#make sure things look right
n = 0
for key, value in alleledict.iteritems():
	print key, value
	n += 1
	if n == 10:
		break

print alleledict['S1_8210']
print len(alleledict.keys())

#since everything looks right, dump out the dict
pickle.dump(alleledict, open("../data/Namfreqs/namfreqdict.pkl", 'wb'))				




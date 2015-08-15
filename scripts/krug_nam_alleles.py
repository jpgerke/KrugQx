
import csv

krugfile = csv.DictReader(open("../data/Krug_freqs_olap.csv", 'rb'))

print krugfile.fieldnames
#build the dict of Krug SNPs
krugdict = {}
for entry in krugfile:
	krugdict[entry['SNP']] = entry

namfile = csv.DictReader(open("../data/NAM_freqs_olap_krug.csv", 'rb'))

#counter of tracked alleles
tracked = 0
allsnps = 0

#files for the output
krugfields = ["SNP","AlleleTracked","A1","A2","KC0Allele1","KLSAllele1","KSSAllele1"]
namwriter =  csv.DictWriter(open("../data/nam_olap_krug_biallelelic.csv", 'wb'), fieldnames = namfile.fieldnames)
krugwriter = csv.DictWriter(open("../data/krug_olap_biallelic.csv", 'wb'), fieldnames = krugfields)
namwriter.writeheader()
krugwriter.writeheader()

#run through the nam file and examine the krug counterpart
for entry in namfile:
	allsnps += 1
	krugsnp = krugdict[entry['SNP']]
	namalleles = [entry['A1'], entry['A2']]
	#is the allele tracked in krug in nam?
	if krugsnp['AlleleTracked'] in namalleles:
		tracked += 1
		#if the alleles are 'backwards' we need to change the frequencies
		if krugsnp['AlleleTracked'] == entry['A2']:
			krugsnp['KC0Allele1'] = 1 - float(krugsnp['KC0Allele1'])
			krugsnp['KLSAllele1'] = 1 - float(krugsnp['KLSAllele1'])
			krugsnp['KSSAllele1'] = 1 - float(krugsnp['KSSAllele1'])
		#unnecessary and redundant sanity check
		else:
			assert krugsnp['AlleleTracked'] == entry['A1']
		#now write out the results
		krugsnp['A1'] = entry['A1']
		krugsnp['A2'] = entry['A2']
		namwriter.writerow(entry)
		krugwriter.writerow(krugsnp)

print tracked
print allsnps

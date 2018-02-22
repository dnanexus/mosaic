#!/usr/env Python
# parse_NCBI_ids.py
# This script processes the abundance table, whereby first column is the organism name

# imports
import sys, os, csv
from ete3 import NCBITaxa			# note: ete3 is best installed with Anaconda/Miniconda
import pandas as p

# updating the taxonomy, will take a couple minutes if running for the first time
ncbi = NCBITaxa()

# setting up definitions and labels
dsets={ 'bio': 'Mouse',
	'sim_low': 'Sim. Low Compl.',
	'sim_medium': 'Sim. Medium Compl.',
	'sim_high': 'Sim. Low Compl.' }
taxonomy=['superkingdom','phylum','class','order','family','genus','species','strain' ]
default={ 'superkingdom': 'k__unknown',
		'phylum': 'p__unknown',
		'class': 'c__unknown',
		'order': 'o__unknown',
		'family': 'f__unknown',
		'genus': 'g__unknown',
		'species': 's__unknown',
		'strain': 'str__unknown' }
abbr={ 'superkingdom': 'k__',
		'phylum': 'p__',
		'class': 'c__',
		'order': 'o__',
		'family': 'f__',
		'genus': 'g__',
		'species': 's__',
		'strain': 'str__' }

dataset=sys.argv[1]
if dataset not in dsets:
	print 'Dataset name "%s" is not valid' % dataset
	exit()
truth_file=sys.argv[2]
subm_file=sys.argv[3]

truth_table=p.read_csv(truth_file, delimiter="\t", header=None)
table=p.read_csv(subm_file, delimiter="\t", header=None)

def parse_table(table,dset_type):
	print "Parsing %s file" % (dset_type)
	outlist=[]
	out=open("profiling_%s_%s_abundances.tsv" % (dataset,dset_type), 'w')
	fileout_names=[]
	for x in range(1,5):
		filename = "%s_%x.tsv" % (dset_type,x)
		outhandle=open(filename,'w')
		fileout_names.append(filename)
		outlist.append(outhandle)

	for i,row in table.iterrows():
		annotation = dict(default)
		tax=row[0].split(",")
		full=[]
		fullname=[]
		for i,id in enumerate(map(int,tax)):
			if i == 7 and int(tax[6]) != 0:
				strain= ncbi.get_taxid_translator([int(tax[6])])[int(tax[6])] + " strain " + str(id)
				full.append("strain")
				fullname.append(strain)
				annotation['strain']=[strain][0]
			elif i == 7 and int(tax[6]) == 0:
				strain= 'str__unknown'
				full.append("strain")
				fullname.append(strain)
				annotation['strain']=[strain][0]
			elif id != 0:
				name=ncbi.get_taxid_translator([id])[id]
				rank=ncbi.get_rank([id])[id]
				full.append(rank)
				fullname.append(name)
				annotation[rank]=[name][0]
		line=[]
		for key in taxonomy:
			if '__unknown' not in annotation[key]:
				line.append(abbr[key]+annotation[key])
			else:
				line.append(annotation[key])

		line.extend(map(str,row)[1:5])
		out.write("\t".join(line)+"\n")
		for n,value in enumerate(row[1:5]):
			if value > 0.0:
				outlist[n].write(str(value)+"\t"+"\t".join(line[0:8])+"\n")
	out.close()
	for outhandle in outlist:
		outhandle.close()
	print "Parsing of %s table complete" % (dset_type)
	return fileout_names

truthfn=parse_table(truth_table,'truth')
subfn=parse_table(table,'submission')

argument=[]
for series in range(1,5):
	idx=series-1
	argument.append('%s,"%s Truth Sample %s"' %  (truthfn[idx],dsets[dataset],series))
	argument.append('%s,"%s Submission Sample %s"' %  (subfn[idx],dsets[dataset],series))
argument_text=" ".join(argument)
cmd="ktImportText %s -o profiling_%s_krona.html" % (argument_text,dataset)
print cmd

os.system(cmd)

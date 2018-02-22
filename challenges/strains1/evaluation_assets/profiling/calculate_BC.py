#!/usr/env Python
# calculate_BC.py
# usage: compare.py one_in two_in col_index
# where one_in and two_in are the abundance tables in terms of the standardized format
# col_index refers to the column number of the sample to be compared in case there are
# multiple samples

import sys, os, scipy
import numpy as np
from scipy import spatial
from collections import defaultdict
from tabulate import tabulate

def compute_jaccard_index(set_1, set_2):
 	n = len(set_1.intersection(set_2))
	return n / float(len(set_1) + len(set_2) - n)

def reformat_infile(infilename,outfilename,col_index):
	infile=open(infilename,'r')
	outfile=open(outfilename,'w')

	for line in infile:
		tmp=line.split("\t")
		tax=['0'] * 7

		index=0
		for i in tmp[0].split(",")[0:-1]:
			tax[index]=i
			index+=1

		taxids=",".join(tax)
		str=",".join([taxids, tmp[col_index]])
		outfile.write(str+"\n")

dataset=sys.argv[1]
one_in=sys.argv[2]
two_in=sys.argv[3]

output=open("profiling_"+dataset+"_braycurtis.tsv", 'wt')
headers=['dataset','sample','tax_ranking','braycurtis']
output.write("\t".join(headers))
output.write("\n")
scores={}
scores["family"]=[]
scores["genus"]=[]
scores["species"]=[]
for col_index in range(0,4):
	print "Processing Sample %s" % str(col_index+1)
	reformat_infile(one_in, one_in+".csv", col_index)
	reformat_infile(two_in, two_in+".csv", col_index)

	one=one_in+".csv"
	two=two_in+".csv"
	clades = [ "kingdom", "phylum", "class", "order" , "family" ,"genus" ,"species"]
	table = []

	for ranking, clade in enumerate(clades):
		taxa = {}
		taxa["one"]= {}
		taxa["two"]= {}
		taxa_n = {}
		taxa_n["one"]= {}
		taxa_n["two"]= {}

		list_of_clades_one = []
		with open(one, 'rt') as f:
			for line in f:
				line = line.rstrip()
				line = line.split(",")
				list_of_clades_one.append(line[ranking])

				if line[ranking] in taxa["one"]:
					taxa["one"][line[ranking]] += float(line[7])
					taxa_n["one"][line[ranking]] += 1
				else:
					taxa["one"][line[ranking]] = float(line[7])
					taxa_n["one"][line[ranking]] = 1
					taxa["two"][line[ranking]] = 0.0
					taxa_n["two"][line[ranking]] = 0

		list_of_clades_two = []
		with open(two, 'rt') as f:
			for line in f:
				line = line.rstrip()
				line = line.split(",")
				list_of_clades_two.append(line[ranking])
				if line[ranking] in taxa["two"]:
					taxa["two"][line[ranking]] += float(line[7])
					taxa_n["two"][line[ranking]] += 1
				else:
					taxa["one"][line[ranking]] = 0.0
					taxa["two"][line[ranking]] = float(line[7])
					taxa_n["one"][line[ranking]] = 0
					taxa_n["two"][line[ranking]] = 1
		list_one = []
		list_two = []
		list_one_n = []
		list_two_n = []
		for key in taxa["one"]:
			list_one.append(taxa["one"][key])
			list_two.append(taxa["two"][key])
			list_one_n.append(taxa_n["one"][key])
			list_two_n.append(taxa_n["two"][key])

		set_one = set(list_of_clades_one)
		set_two = set(list_of_clades_two)
		jaccard= compute_jaccard_index(set_one,set_two)

		sim = 1 - scipy.spatial.distance.braycurtis(list_one,list_two)
		sim_n = 1 - scipy.spatial.distance.braycurtis(list_one_n,list_two_n)
		table.append([clade,sim,sim_n,sim*sim_n,jaccard])

	print tabulate(table, headers=["Clade","Bray Curtis Abundance", "Bray Curtis OTU", "BC Ab * BC OTU", "Jaccard Similarity"], tablefmt="fancy_grid")

	for i in range(4,7):
		scores[table[i][0]].append(table[i][1])
		print [dataset, str(col_index+1),table[i][0],table[i][1]]
		output.write("\t".join([dataset, str(col_index+1), table[i][0], str(table[i][1])] ))
		output.write("\n")

	os.remove(one)
	os.remove(two)
output.close()

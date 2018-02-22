#!/usr/env Python
# compare_results.py
# Created by Sam Westreich, swestreich@dnanexus.com, github.com/transcript
# This script now sets defined thresholds for cutoff values, to increase overall speed and better allow comparisons between different submission results.

# imports
import sys, math
import numpy as np
from sklearn import metrics as skmetrics

# starting files
print "profiling input type should be ARGV1 (sim_low, sim_med, sim_high, or biological), truth file should be ARGV2, submission file should be ARGV3."
truth_file = sys.argv[2]
results_file = sys.argv[3]

def read_answer_key(file):
	truth_dictionary = {}
	with open(file, 'r') as infile:
		for line in infile:
			line = line.strip().split("\t", 1)
			truth_dictionary[line[0]] = map(float, line[1].split("\t"))
	return truth_dictionary

def read_submission(file):
	submission_dictionary = {}
	line_count = 0
	with open(file, 'r') as infile:
		for line in infile:
			line_count += 1
			line = line.strip().split("\t", 1)
			entry_matrix = map(float, line[1].split("\t"))
			# some sanity checking
			if len(entry_matrix) != 4:
				sys.exit("Warning: wrong number of columns in submission file.")
			submission_dictionary[line[0]] = entry_matrix
	return submission_dictionary

def get_stats_strain(truth_dictionary, submission_dictionary, cutoff):
	submission_strains = submission_dictionary.keys()
	for entry in submission_strains:
		if entry.rsplit(',',1)[1] == "0":
			submission_strains.remove(entry)
	truth_strain = truth_dictionary.keys()
	for entry in truth_strain:
		if entry.rsplit(',',1)[1] == "0":
			truth_strain.remove(entry)
	TP = 0
	removed_count = 0
	for tax_entry in submission_strains:
		# thresholding
		cut = True
		for value in submission_dictionary[tax_entry]:
			if value > cutoff:
				cut = False
		if cut == True:
			removed_count += 1
			continue						# skipping rows below cutoff
		# if we've made it here, we're above the cutoff threshold
		if tax_entry in truth_strain:
			TP += 1
	FP = len(submission_strains) - TP - removed_count 	# every strain called incorrectly in submission
	FN = len(truth_strain) - TP 						# every strain missed in submission
	return TP, FP, FN

def get_stats_species(truth_dictionary, submission_dictionary, cutoff):
	# thresholding
	filtered_list = []
	removed_count = 0
	truth_removed_count = 0
	for tax_entry in submission_dictionary.keys():
		cut = True
		for value in submission_dictionary[tax_entry]:
			if value > cutoff:
				cut = False
		if cut == True:
			removed_count += 1
			filtered_list.append(tax_entry.rsplit(',',1)[0])	# this removes the last entry (strain)
		# removing entries with no species info
		elif tax_entry.rsplit(',',1)[1] == "0":
			removed_count += 1
			filtered_list.append(tax_entry.rsplit(',',1)[0])
	# now with the thresholded entries:
	submission_species = set(entry.rsplit(',',1)[0] for entry in submission_dictionary.keys())
	truth_species = set(entry.rsplit(',',1)[0] for entry in truth_dictionary.keys())
	for entry in set(truth_species):
		if entry.rsplit(',',1)[1] == "0":
			truth_removed_count += 1
	TP = 0
	for tax_entry in submission_species:
		if tax_entry not in filtered_list and tax_entry in truth_species:
			TP += 1
	FP = len(submission_species) - TP - removed_count 		# every strain called incorrectly in submission
	if FP < 0:
		FP = 0
	FN = len(truth_species) - TP - truth_removed_count			# every strain missed in submission
	return TP, FP, FN

def get_stats_genus(truth_dictionary, submission_dictionary, cutoff):
	# thresholding
	filtered_list = []
	removed_count = 0
	for tax_entry in submission_dictionary.keys():
		cut = True
		for value in submission_dictionary[tax_entry]:
			if value > cutoff:
				cut = False
		if cut == True:
			removed_count += 1
			filtered_list.append(tax_entry.rsplit(',',2)[0])	# removes the last two entries (species & strain)
	# now with the thresholded entries:
	submission_genus = set(entry.rsplit(',',2)[0] for entry in submission_dictionary.keys())
	truth_genus = set(entry.rsplit(',',2)[0] for entry in truth_dictionary.keys())
	TP = 0
	for tax_entry in submission_genus:
		if tax_entry not in filtered_list and tax_entry in truth_genus:
			TP += 1
	FP = len(submission_genus) - TP - removed_count 		# every strain called incorrectly in submission
	if FP < 0:
		FP = 0
	FN = len(truth_genus) - TP 			# every strain missed in submission
	return TP, FP, FN

def compute_metrics(metrics_list):
	TP = float(metrics_list[0])			# true positives
	FP = float(metrics_list[1])			# false positive - in submission but NOT in truth file
	FN = float(metrics_list[2])			# false negatives - in truth file but NOT in submission
	try:
		precision = TP/(TP+FP)
	except ZeroDivisionError:
		precision = 1.0
	recall = TP/(TP+FN)
	try:
		F1_score=2*(precision*recall)/(precision+recall)
	except ZeroDivisionError:
		F1_score = 0
	misclass=(FP+FN)/(TP+FP+FN)
	return precision, recall, F1_score

def iterate_loop(submission_dictionary, truth_dictionary, level, iter_outfile):
	# set up all the different confidence thresholds
	iterate_values = []
	iterate_starting_values = [0.000001, 0.000002, 0.000003, 0.000004, 0.000005, 0.000006, 0.000007, 0.000008, 0.000009]
	for val in iterate_starting_values:
		iterate_values.append(val)
		iterate_values.append(val*10.0)
		iterate_values.append(val*100.0)
		iterate_values.append(val*1000.0)
		iterate_values.append(val*10000.0)
		iterate_values.append(val*100000.0)
	iterate_values = sorted(iterate_values, key=float)
	precision_list = []
	recall_list = []
	max_f1_score = 0							# used for tracking the highest F1 score and cutoff to get that score
	max_f1_cutoff = 0
	iter_outfile.write("cutoff\tTP\tFN\tFP\tPrecision\tRecall\tF1\n")
	for val in iterate_values:
		if level == "strain":
			stats = get_stats_strain(truth_dictionary, submission_dictionary, val)
			if stats[0] == 0:			# no more true positives
				break
			iter_outfile.write(str(val) + "\t" + str(stats[0]) + "\t" + str(stats[2]) + "\t" + str(stats[1]) + "\t")
			metrics = compute_metrics(stats)
			iter_outfile.write(str(metrics[0]) + "\t" + str(metrics[1]) + "\t" + str(metrics[2]) + "\n")
		elif level == "species":
			stats = get_stats_species(truth_dictionary, submission_dictionary, val)
			if stats[0] == 0:			# no more true positives
				break
			iter_outfile.write(str(val) + "\t" + str(stats[0]) + "\t" + str(stats[2]) + "\t" + str(stats[1]) + "\t")
			metrics = compute_metrics(stats)
			iter_outfile.write(str(metrics[0]) + "\t" + str(metrics[1]) + "\t" + str(metrics[2]) + "\n")
		elif level == "genus":
			stats = get_stats_genus(truth_dictionary, submission_dictionary, val)
			if stats[0] == 0:			# no more true positives
				break
			iter_outfile.write(str(val) + "\t" + str(stats[0]) + "\t" + str(stats[2]) + "\t" + str(stats[1]) + "\t")
			metrics = compute_metrics(stats)
			iter_outfile.write(str(metrics[0]) + "\t" + str(metrics[1]) + "\t" + str(metrics[2]) + "\n")
		else:
			print "Level not properly provided."
		if metrics[0] > 0.0:
			precision_list.append(metrics[0])
			recall_list.append(metrics[1])
		# for calculating max F1 score, we keep track through the different cutoffs to see which cutoff gives the highest score.
		if metrics[2] > max_f1_score:
			max_f1_score = metrics[2]
			max_f1_cutoff = val
	auprc = skmetrics.auc(recall_list, precision_list)		# area under precision/recall curve
	return max_f1_score, max_f1_cutoff, auprc

# reading in the input files
truth_dic = read_answer_key(truth_file)
submission_dic = read_submission(results_file)
dataset = sys.argv[1]

# creating the precision-recall curve results for STRAIN at each cutoff threshold
strains_outfile = open("profiling_" + dataset + "_PRC_strain.tsv", 'w')
strains_iteration = iterate_loop(submission_dic, truth_dic, "strain", strains_outfile)
strains_outfile.close()

# creating the precision-recall curve results for SPECIES at each cutoff threshold
species_outfile = open("profiling_" + dataset + "_PRC_species.tsv", 'w')
species_iteration = iterate_loop(submission_dic, truth_dic, "species", species_outfile)
species_outfile.close()

# creating the precision-recall curve results for GENUS at each cutoff threshold
genus_outfile = open("profiling_" + dataset + "_PRC_genus.tsv", 'w')
genus_iteration = iterate_loop(submission_dic, truth_dic, "genus", genus_outfile)
genus_outfile.close()

# writing the final scores outfile
score_outfile = open("profiling_" + dataset + "_scores.tsv", "w")
headers = ["tax_ranking", "dataset", "TP", "FN", "FP", "Precision", "Recall", "F1", "improved_F1", "cutoff", "AUC"]
score_outfile.write("\t".join(headers) + "\nstrain\t" + dataset + "\t")
# writing the row for strains...
strain_stats = get_stats_strain(truth_dic, submission_dic, 0)
score_outfile.write(str(strain_stats[0])+"\t"+str(strain_stats[2])+"\t"+str(strain_stats[1]) + "\t")				# TP, FN, FP
print strain_stats
strain_metrics = compute_metrics(strain_stats)
score_outfile.write("\t".join(str(x) for x in strain_metrics) + "\t")			# precision, recall, F1
score_outfile.write("\t".join(str(x) for x in strains_iteration) + "\n")		# improved_F1, cutoff, AUC
# ...for species...
score_outfile.write("species\t" + dataset + "\t")
species_stats = get_stats_species(truth_dic, submission_dic, 0)
score_outfile.write(str(species_stats[0])+"\t"+str(species_stats[2])+"\t"+str(species_stats[1]) + "\t")				# TP, FN, FP
print species_stats
species_metrics = compute_metrics(species_stats)
score_outfile.write("\t".join(str(x) for x in species_metrics) + "\t")			# precision, recall, F1
score_outfile.write("\t".join(str(x) for x in species_iteration) + "\n")		# improved_F1, cutoff, AUC
# ...and for genus...
score_outfile.write("genus\t" + dataset + "\t")
genus_stats = get_stats_genus(truth_dic, submission_dic, 0)
score_outfile.write(str(genus_stats[0])+"\t"+str(genus_stats[2])+"\t"+str(genus_stats[1]) + "\t")				# TP, FN, FP
print genus_stats
genus_metrics = compute_metrics(genus_stats)
score_outfile.write("\t".join(str(x) for x in genus_metrics) + "\t")			# precision, recall, F1
score_outfile.write("\t".join(str(x) for x in genus_iteration) + "\n")		# improved_F1, cutoff, AUC
score_outfile.close()

print "Run successful."		# success!

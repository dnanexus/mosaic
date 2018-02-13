#!/usr/bin/env python

# matrix_strains2_evaluator.py
# Created 1/26/18 by Sam Westreich

# Purpose:
#	1. Read in truth file, submission file
#	2. Compare and calculate accuracy, precision, F1, misclass, ARI
#	3. Calculate matrix distance
#	4. Loop of:
#		1. Remove lowest confidence prediction
#		2. Recalculate accuracy, precision, F1, misclass, ARI
#		3. Save ARI and removed confidence level if ARI > maximum
#		4. Go to 1 (x39)
#	5. Return #2 as file output 1, return #4 results as file output 2

# Imports
import sys, os, math
import sklearn
from sklearn import metrics
from sklearn.metrics import f1_score
from sklearn.metrics.cluster import adjusted_rand_score
import numpy as np

# starting files
truth_file=sys.argv[1]
results_file=sys.argv[2]

def read_answer_key(file):
	answer_matrix = []
	strains_list = []
	with open(file, 'r') as infile:
		for line in infile:
			line = line.strip().split("\t", 1)
			strains_list.append(line[0])
			matrix_entry = map(float, line[1].split("\t"))

			answer_matrix.append(matrix_entry)

	return answer_matrix, strains_list

def read_submission(file):
	submission_matrix = []
	line_count = 0
	with open(file, 'r') as infile:
		for line in infile:
			line_count += 1
			line = line.strip().split("\t", 1)
			entry_matrix = map(float, line[1].split("\t"))

			# some sanity checking
			if len(entry_matrix) != 4:
				sys.exit("Warning: wrong number of columns in submission file.")
			test_sum = 0
			for entry in entry_matrix:
				test_sum += math.ceil(float(entry))
			if test_sum > 1.0:
				sys.exit("Warning: organism marked as present in more than 1 sample.")

			submission_matrix.append(entry_matrix)

	# another sanity check
	if line_count != 40:
		sys.exit("Warning: incorrect number of lines in submission file.")

	return submission_matrix

def get_stats(truth_matrix, submission_matrix):
	TP = 0.0	# there's a number where the truth is nonzero
	FP = 0.0	# there's a number where the truth is zero
	FN = 0.0	# there's no number where the truth is nonzero
	TN = 0.0	# there's no number where the truth is zero

	# we multiply the answer key by 10, to avoid some subtraction issues
	truth_matrix = np.multiply(truth_matrix, 10.0)

	# we subtract the answers from the answer key
	difference_matrix = truth_matrix - submission_matrix

	# we can now evaluate for TP, FP, FN, TN
	for row in difference_matrix:
		if sum(row) == 0.0:
			# It was all 0s in both the answer key and the submission
			TN += 1
		elif sum(row) == 10.0:
			# Entry in answer key but not in submission
			FN += 1
		elif 0.0 < sum(row) < 10.0:
			if np.count_nonzero(row) == 1:
				# Three 0s, one entry of less than ten means that the answer is right!
				TP += 1
			else:
				# They got the wrong sample for this strain
				FP += 1
		else:
			# row sum is less than 0, the submission thinks a strain is present where it isn't
			FP += 1

	return TP, FP, TN, FN

def compute_metrics(metrics_list):
	# this part is totally taken from Keng's data, thanks Keng
	TP = metrics_list[0]
	FP = metrics_list[1]
	TN = metrics_list[2]
	FN = metrics_list[3]
	total=TP+TN+FP+TN
	accuracy=(TP+TN)/total
	precision=TP/(TP+FP)
	recall=TP/(TP+FN)
	F1_score=2*(precision*recall)/(precision+recall)
	misclass=(FP+FN)/total								# fixed, Keng had this wrong

	return [accuracy,precision,recall,F1_score,misclass]

def adjusted_rand(answers_matrix, submission_matrix):
	# the adjusted_rand_score option just takes a 1D array, so we've got to flatten these 2D arrays.
	answers_list = []
	submission_list = []
	for sublist in answers_matrix:
		for entry in sublist:
			answers_list.append(entry)
	for sublist in submission_matrix:
		for entry in sublist:
			submission_list.append(entry)

	# now we calculate rand score
	return adjusted_rand_score(answers_list, submission_list)

truth_matrix = read_answer_key(truth_file)[0]
submission_matrix = read_submission(results_file)

# creating the first output file
stats_outfile = open(results_file[:-4] + "_scores.tsv", "w")

# header
stats_outfile.write("TP\tFP\tTN\tFN\taccuracy\tprecision\trecall\tF1_score\tmisclassification_rate\tadjusted_rand_index\n")
# data
init_metrics = compute_metrics(get_stats(truth_matrix, submission_matrix))
stats_outfile.write("\t".join(str(int(item)) for item in get_stats(truth_matrix, submission_matrix)) + "\t")
stats_outfile.write("\t".join(str(item) for item in init_metrics))
stats_outfile.write("\t" + str(adjusted_rand(truth_matrix, submission_matrix)))
stats_outfile.close()

# Now, we need to generate the second output file...

# sanity check - is this file binary?
flat_submission_list = [item for sublist in submission_matrix for item in sublist]
if list(set(flat_submission_list)) == [0,1]:
	binary_report = open(results_file[:-4] + "_binary", "w")
	binary_report.write("binary == true")
	binary_report.close()
	sys.exit()

# This one will require an iteration of removing the lowest confidence score, over and over.
iter_outfile = open(results_file[:-4] + "_PRC.tsv", "w")
iter_outfile.write("cutoff\tTP\tFP\tTN\tFN\taccuracy\tprecision\trecall\tf1\tmisclassification\tadj_rand_index\n")

for iteration in range(0,39):
	# first, we get the lowest number in the submission
	masked_submission_matrix = submission_matrix[submission_matrix != 0.0]
	lowest = np.min(masked_submission_matrix)

	# second, we remove the row with this lowest confidence value
	row_counter = 0
	for row in submission_matrix:
		if lowest in row:
			submission_matrix = np.delete(submission_matrix, row_counter, axis=0)
			truth_matrix = np.delete(truth_matrix, row_counter, axis=0)
			if lowest == 0.0:
				break
		else:
			# we need to continue in case there are multiple instances of the same cutoff
			row_counter += 1
	# At this point, if the lowest value is 0 confidence, we skip it (go to next loop iteration)
	if lowest == 0.0:
		continue
	# third, we recalculate our stats
	try:
		stats = get_stats(truth_matrix, submission_matrix)
		metrics = compute_metrics(get_stats(truth_matrix, submission_matrix))
		iter_outfile.write(str(lowest) + "\t" + "\t".join(str(int(item)) for item in stats) + "\t")
		iter_outfile.write("\t".join(str(item) for item in metrics) + "\t" + str(adjusted_rand(truth_matrix, submission_matrix)) + "\n")
	except ZeroDivisionError:
		# this occurs when trying to divide by 0, obviously
		# At this point, we're out of positive values to subtract.
		break

iter_outfile.close()

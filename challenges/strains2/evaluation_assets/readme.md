# Strains #2 challenge evaluation script

***

**Contents:**

1. strains2_evaluator.py

***

**Python package dependencies:**

1. sys
2. os
2. math
3. numpy
4. sklearn (f1_score, adjusted_rand_score)

***

**Usage:**

    $ python strains2_evaluator.py truth_file.txt submission_file.txt
    
***

**How it works:**

This Python scripts evaluates submissions to the Strains #2 challenge, grading them against the truth file.  It returns two files: the score of the submission file as presented and, if the file contains confidence scores (0<score<1 for each entry), it calculates the scores at each confidence threshold level.  This allows for construction of the precision/recall curve that's presented on the Mosaic challenge site.

The script performs the following steps:

1. Read in the answers from the truth file.
2. Reads in the submission file.
3. Compares the two files to determine the number of true positives, false positives, true negatives, and false negatives.
	1. True positive = strain is present in both submission and answer key in the proper metagenome sample.
	2. False positive = strain marked as present in the submission file in a sample where it is not present, according to the truth file.
	3. True negative = strain is correctly marked as not present for all four metagenomes in the submission file.
	4. False negative = strain is marked as not present for all four metagenomes in the submission file, when it is, in fact, present in one of the metagenomes, according to the truth file.
4. Uses the true positives, false positives, true negatives, and false negatives to calculate accuracy, precision, recall, and F1 score.
5. Compares the two files to create the adjusted Rand index score.
6. If the submission file contains confidence estimates, iterates over the submission file, dropping the lowest confidence estimate and rescoring.  This continues until the file has been scored at each confidence threshold.
7. If the submitted file is binary, an empty outfile named $filename_binary.txt is created as output instead of step 6.

***

**Outputs:**

$submission-filename_scores.tsv - contains scores for submission file.

$submission-filename_PRC.tsv - contains the scores at each confidence cutoff, used in building the precision/recall curve.

OR

$submission-filename_binary - empty file, marks the submission as being binary (no confidence scores given).
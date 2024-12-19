#!/usr/bin/env python3

# parse_TRF.1.py

# Shu-Ting Cho <vivianlily6@hotmail.com>
# Read TRF's result and output tsv

# v1 2023/10/25

# Usage:
# python3 parse_TRF.1.py --in_file=/Users/shc167/Documents/project/TR02/TR02.89/KEGG_nts_v2.trf/K23578.tsv --out_file=/Users/shc167/Documents/project/TR02/TR02.89/TRF_tsv/K23578.tsv

## import modules
#import re
import argparse

# parsing arguments
def parse_args():
	parser = argparse.ArgumentParser(description="Read predict result and output tsv.")
	parser.add_argument('-i', '--in_file', help="predict result")
	parser.add_argument('-o', '--out_file', help="output .tsv file")
	return parser.parse_args()

args = parse_args()

out_f = open(str(args.out_file), "w")
out_f.write("seqID\tstart\tend\tPeriodSize\tCopyNumber\tConsensusSize\tPercentMatches\tPercentIndels\tScore\tA\tC\tG\tT\tEntropy\n")

# Opening result file
seqIDs = []
with open(args.in_file) as in_f:
	lines = in_f.readlines()
	for line in lines:
		line = line.rstrip()
		if line[0] == '@':
			seqID = line[1:]
		elif line[0].isdigit():
			out_words = [seqID] + line.split()[0:13]
			out_f.write('\t'.join(out_words) + '\n')
			seqIDs.append(seqID)
out_f.close()

print('n_seq_in =', len(set(seqIDs)))
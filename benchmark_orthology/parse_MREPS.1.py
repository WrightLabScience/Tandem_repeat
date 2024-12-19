#!/usr/bin/env python3

# parse_MREPS.1.py

# Shu-Ting Cho <vivianlily6@hotmail.com>
# Read MREPS's result and output tsv

# v1 2023/10/25

# Usage:
# python3 parse_MREPS.1.py --in_file=/Users/shc167/Documents/project/TR02/TR02.89/KEGG_nts_v2.MREPS/K23578.tsv --out_file=/Users/shc167/Documents/project/TR02/TR02.89/MREPS_tsv/K23578.tsv

## import modules
import re
import argparse

# parsing arguments
def parse_args():
	parser = argparse.ArgumentParser(description="Read predict result and output tsv.")
	parser.add_argument('-i', '--in_file', help="predict result")
	parser.add_argument('-o', '--out_file', help="output .tsv file")
	return parser.parse_args()

args = parse_args()

out_f = open(str(args.out_file), "w")
out_f.write("seqID\tstart\tend\tsize\tperiod\texponent\terr\n")

# Opening result file
pattern_hit = re.compile('^\s+[0-9]+.*$') # starts from spaces and numbers

seqIDs = []
with open(args.in_file) as in_f:
	lines = in_f.readlines()
	for line in lines:
		line = line.rstrip()
		if line[0:19] == 'Processing sequence':
			seqID = line[21:101]
		elif pattern_hit.match(line):
			hit_words = line.split()
			start, end, size, err = hit_words[0], hit_words[2], hit_words[4], hit_words[7]
			period = hit_words[5].strip('<>')
			exponent = hit_words[6].strip('[]')
			out_f.write('\t'.join([seqID, start, end, size, period, exponent, err]) + '\n')
			seqIDs.append(seqID)
out_f.close()

print('n_seq_in =', len(set(seqIDs)))
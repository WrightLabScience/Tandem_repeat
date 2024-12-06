#!/usr/bin/env python3

# parse_RADAR.1.py

# Shu-Ting Cho <vivianlily6@hotmail.com>
# Read RADAR's result and output tsv

# v1 2023/11/16

# Usage:
# python3 parse_RADAR.1.py --in_dir=/Users/shc167/Documents/project/TR02/TR02.93/RADAR.TPout/ --out_file=/Users/shc167/Documents/project/TR02/TR02.93/RADAR.TPout.tsv

## import modules
import re, os
import argparse

# parsing arguments
def parse_args():
	parser = argparse.ArgumentParser(description="Read predict result and output tsv.")
	parser.add_argument('-i', '--in_dir', help="predict result")
	parser.add_argument('-o', '--out_file', help="output .tsv file")
	return parser.parse_args()

args = parse_args()

out_f = open(str(args.out_file), "w")
out_f.write("seqID\tstart\tend\tscore\tNo_of_Repeats\tTotal_Score\tLength\tDiagonal\tBW_From\tBW_To\tLevel\n")

# Opening result file
pattern_hit = re.compile('^\s+([0-9]+)-\s+([0-9]+)\s.*$') # starts from spaces and numbers

seqIDs = []

for file in os.listdir(args.in_dir):
	# check only output files
	if file.endswith('.out.txt'):
		seqID = file.split('.')[0]
		full_path = args.in_dir.rstrip('/') + '/' + file
		with open(full_path) as in_f:
			lines = in_f.readlines()
			count_area = 0
			for line in lines:
				line = line.rstrip()
				if line.startswith('No repeats found'):
					continue
				elif count_area<2 and line.startswith('---'):
					count_area += 1
				elif count_area==1 and not line.startswith('No. of Repeats'):
					No_of_Repeats, Total_Score, Length, Diagonal, BW_From, BW_To, Level = line.replace(" ", "").split('|')
					start, end = 'start', 0
				elif count_area==2 and not line.startswith('---'):
					range_h = pattern_hit.search(line)
					start_new, end_new = range_h.group(1), range_h.group(2)
					if start=='start':
						start = start_new
					if int(end_new) > int(end):
						end = end_new
				elif count_area==2 and line.startswith('---'):
					count_area = 0 # end of the record
					out_f.write('\t'.join([seqID, start, end, Total_Score, No_of_Repeats, Total_Score, Length, Diagonal, BW_From, BW_To, Level]) + '\n')
		seqIDs.append(seqID)
out_f.close()

print('n_seq_in =', len(set(seqIDs)))
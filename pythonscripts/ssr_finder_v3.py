#!/usr/bin/env python

from __future__ import print_function
from Bio import SeqIO
from Bio.Seq import Seq
import regex as re
import argparse
from tqdm import tqdm
import sys

# User input functions
parser = argparse.ArgumentParser(prog="SSR finder program")
parser.add_argument("-fi", dest = "fa", required=True, help = "Fasta file in which the repeats have to be located")
parser.add_argument("-rep", dest = "rep", required=True, help = "Repeat file which contains all the repeat sequences")
parser.add_argument("-min", "--repeat_len", type = int, dest = "min", default = 12, help = "The minimum repeat length to be considered")
# parser.add_argument("-mut", "--mutation", type = int, dest = "mutation", help = "The minimum number of mutations in the repeat to be considered")
# parser.add_argument("-sub", "--substitution", type = int, dest = "mutation", help = "The minimum number of mutations in the repeat to be considered")
# parser.add_argument("-indel", "--mutation", type = int, dest = "mutation", help = "The minimum number of mutations in the repeat to be considered")
parser.add_argument("-fo","--out", type = argparse.FileType('w'), default = sys.stdout, help = "The output file in which the results will be written")
args = parser.parse_args()

input_file = args.fa
repeat_file = args.rep
output_file = args.out
min_length = args.min


# Function for building a regular expression
def build_regex(pattern):
	suffix = []
	prefix = []
	length = len(pattern)
	if length == 1:
		regex = "("+pattern+"){12,}"
	else:
		n = int((12/length))-1
		for i in reversed(range(1,length)):
			prefix.append(pattern[-i:])
			suffix.append(pattern[:i])
		suffix = "("+("|".join(suffix))+"){0,1}"
		prefix = "("+("|".join(prefix))+"){0,1}"
		center = "("+pattern+"){"+str(int(n))+",}"
		regex = prefix + center + suffix
	return regex

#List of repeats from the file
repeats = []
with open(repeat_file) as rep_file:
	for line in rep_file:
		rep = line.strip()
		repeats.append(rep)	

seqs = list(SeqIO.parse(input_file, "fasta"));

for record in seqs:
	seq_id = record.id
	seq = (str(record.seq)).upper()
	seq_len = len(seq)
	for repeat in tqdm(repeats ,total=len(repeats), desc = "Processing " + seq_id + " " + str(seq_len) + " bp ", leave = True, unit = "Repeats", position = 1):
		rclass = repeat
		rlen = len(repeat)
		
		revcom = str(Seq(repeat).reverse_complement())
		if revcom == rclass:
			patterns = [rclass]
		elif revcom != rclass:
			patterns = [rclass, revcom]

		for p, pattern in enumerate(patterns):
			regex = re.compile(build_regex(pattern))
			
			if p == 0:
				strand = '+'
			else:
				strand = '-'
			
			for a in regex.finditer(seq):
				ssr_start = a.start()
				ssr_stop = a.end()
				ssr_seq = a.group(0) 
				ssr_len = len(ssr_seq)
				up_seq = seq[ssr_start-rlen : ssr_start]
				prefix = ssr_seq.find(pattern)
				if prefix == 0:
					missp = pattern
				elif prefix>0:
					missp = pattern[:-prefix]
				for i in range(1, len(missp) + 1):
					j = -i
					if (len(up_seq) == rlen) and (up_seq[j] == missp[j]):
						ssr_start = ssr_start - 1
						ssr_len = ssr_len + 1
						ssr_seq = missp[j] + ssr_seq
					else:
						break
				if ssr_len>= min_length:
					print(record.id, ssr_start + 1, ssr_stop, ssr_len, ssr_seq[:rlen], rclass, sep="\t", file = output_file)

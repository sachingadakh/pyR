#!usr/bin/python

from __future__ import print_function,division
import sys
import pysam
import argparse
import numpy as np
from tqdm import tqdm
from collections import OrderedDict

parser = argparse.ArgumentParser(prog="Bam Coverage Program")
parser.add_argument("bed", help="Bed file for which the coverage should be reported")
parser.add_argument("-bam", dest="bam", nargs='+', required=True, help="Bam files across which the coverage is calculated")
parser.add_argument("-f", "--flank", type=int, dest="flank", default=3000, help="")
parser.add_argument("-w", "--window", type=int, dest="window", default=100, help="")
parser.add_argument("-o", "--out", type=argparse.FileType('w'), default=sys.stdout, help="")
parser.add_argument("-n", "--norm", action='store_true', default=False)
parser.add_argument("-p", "--plot", action='store_true', default=False, help="Plot averages of each file")
args = parser.parse_args()

#Gives the values of windows
def getwindows(flank,window):
	windows = []
	tail = flank % window
	nofwindws = flank // window
	for i in xrange(nofwindws+1):
		windows.append(i*window*(-1))
	for i in xrange(1,nofwindws+1):
		windows.append(i*window)
	if tail != 0:
		windows.insert(nofwindws+1,windows[nofwindws]-tail)
		windows.append(windows[len(windows)-1]+tail)
	else:
		pass
	return windows


windows = getwindows(3000,25)

#Gives the array for division of to get normalised coverages. Output is a numpy array.i
def getnormaliser(bamfiles):
	normaliser = []
	for bamfile in bamfiles:
		samfile = pysam.AlignmentFile(bamfile)
		normaliser.append(int(samfile.mapped))
	min_reads = min(normaliser)
	normaliser = np.array(normaliser)
	if min_reads > 20000000:
		normal = 20000000
	elif min_reads > 10000000:
		normal = 10000000
	elif min_reads > 5000000:
		normal = 5000000
	else:
		normal = 1000000
	normaliser = normaliser/normal
	return normaliser


flank = args.flank
window = args.window
if window == 0:
	windows = [0]
else:
	windows = getwindows(flank,window)
bedfile = args.bed
bamfiles = args.bam
norm = args.norm
plot = args.plot

with open(bedfile) as bedf:
	lines = len(bedf.readlines())

if norm:
	normal = getnormaliser(bamfiles)

coverage_dict = OrderedDict()


for k,bamfile in tqdm(enumerate(bamfiles), total = len(bamfiles)):
	print("\t Processing ",bamfile, file=sys.stderr)
	fname = bamfile[:-4]
	coverages = []
	samfile = pysam.AlignmentFile(bamfile, 'rb')
	with open(bedfile) as bedf:
		for entry in tqdm(bedf,total = lines):
			entry = entry.strip().split()
			ref = entry[0]
			ori_start = int(entry[1])
			ori_stop = int(entry[2])
			coverage = []
			for win in windows:
				if win == 0:
					start = ori_start
					stop = ori_stop
				elif win < 0:
					start = ori_start + win
					stop = start + window
				elif win > 0:
					stop = ori_stop + win
					start = stop - window
				#try:
				cov = samfile.count(reference = ref, start = start, end = stop)
				#except:
					#cov = 0
				if norm:
					cov = round(cov/normal[k],3)
				else:
					pass
				coverage.append(cov)
			if plot:
				coverages.append(coverage)
			entry_new = map(str, entry + [fname] + coverage)
			print('\t'.join(entry_new), file = args.out)
	print("\t done!", bamfile, file=sys.stderr)

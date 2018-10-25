#!/usr/bin/env python

"""
Written for python 2.7 or later. Needs pysam package

Report the coverage of reads for each BED entry. Reads that partially overlap the BED entry are also counted. By default, the coverage is normalized to total number of mapped reads. The scaling factor is dependent on the file with least sequencing depth (read below)

CAUTION: Please use sorted BED input only. As the script works by seeking the bam file for each entry, using a random (unsorted) BED can severely slow it down. In case you MUST use an unsorted BED file, bedtools multicov will be a faster approach

Scaling Factor: The factor is decided based on the bam file with least sequencing depth. Here are the four possibilities:
        Least Sequencing Depth     Scale
        >20 million                20mil
        >10 million                10mil
        >5 million                 5mil
        <5 million                 1mil
If you are sure of the sequencing depth of your files, and want to use a different scale option, set it using the -s/--scale option
"""

from __future__ import print_function, division
import sys
import argparse
from collections import defaultdict
from itertools import takewhile, repeat
import numpy as np
import pysam
from datetime import datetime
from tqdm import tqdm

__version__ = 'v0.2'
# Get command line arguments
def getArgs():
    parser = argparse.ArgumentParser(prog='bam_coverage', formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
    parser.add_argument('-bams', nargs='+', metavar='<FILE>', help='input BAM files. Only sorted and indexed BAM files are accepted')
    parser.add_argument('-bed', required=True, metavar='<FILE>', help='BED file for which coverage should be reported (sorted BED recommended)')
    parser.add_argument('-o', '--out', type=argparse.FileType('w'), default=sys.stdout, metavar='<FILE>', help='output file (default is STDOUT)')
    parser.add_argument('-e', '--err', type=argparse.FileType('w'), default=sys.stderr, metavar='<FILE>', help='report progress to this file (default is STDERR)')
    parser.add_argument('--no-header', action='store_true', default=False, help='do not print header in the output file')
    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)
    #parser.add_argument('--strict', default=False, action='store_true', help='validate BED file before processing, and exit the program if any errors')

    progress = parser.add_mutually_exclusive_group()
    progress.add_argument('-q', '--quiet', action='store_true', default=False, help='do not print any progress')
    progress.add_argument('-v', '--verbose', action='store_true', default=False, help='print progress in more detail')

    norm = parser.add_mutually_exclusive_group()
    norm.add_argument('--no-norm', default=False, action='store_true', help='do not normalize the coverage to sequencing depth')
    norm.add_argument('-s', '--scale', default=None, metavar='INT', type=int, help='scale to be used for normalization (Read program description)')

    return parser.parse_args()

# Returns coverage of a BED entry for a given BAM file
def getCoverage(bamfile, bedEntry):
    l = bedEntry.strip().split('\t')
    ref, start, stop = l[0:3]
    try:
        coverage = bamfile.count(reference=ref, start=int(start), end=int(stop))
    except ValueError:
        coverage = 0
    return coverage

# Returns an array of coverage as normalized to depth. Uses numpy arrays for vectorized division
def normalizeCoverage(array, depth):
    assert len(array[0]) == len(depth), 'Arrays of coverage and depth are not of equal size!'
    return array/depth

# Returns appropriate scaling factor based on sequencing depth of BAM files
def getScale(depthdict):
    minimumDepth = min(depthdict.values())
    if args.verbose:
        print("Minimum sequencing depth in BAM files is %d" % minimumDepth, file=args.err)
    if minimumDepth > 20000000:
        return 20000000
    elif minimumDepth > 10000000:
        return 10000000
    elif minimumDepth > 5000000:
        return 5000000
    else:
        return 1000000

# Self explanatory
def printHeader(outfile, bedfile):
    header = ['#Chr', 'Start', 'Stop', 'Name', 'Score', 'Strand']
    numColumns = 0
    with open(bedfile) as bedFile:
        first_row = bedFile.next()
        numColumns = len(first_row.strip().split('\t'))
    if numColumns > 6:
        for i in range(6, numColumns):
            header.append('ExtraColumn_' + str(i-5))
    header = header[:numColumns] + fileNames
    print(*header, sep="\t", file=outfile)

# Returns line count of a given file. Uses buffered raw reading for very fast computation
def rawLineCount(filename):
    f = open(filename, 'rb')
    bufgen = takewhile(lambda x: x, (f.read(1024*1024) for _ in repeat(None)))
    return sum( buf.count(b'\n') for buf in bufgen if buf )

# Self explanatory
def processBAMFile(file, num=0):
    bamfile = pysam.AlignmentFile(file, 'rb')
    fileNames.append(file)
    mappedReads[file] = bamfile.mapped
    if args.verbose:
        checkPoint = totalEntries//1000
        print("Total number of mapped reads in %s: %d" % (file, bamfile.mapped), file=args.err)
    with open(args.bed, 'r') as bedFile:
        for i, entry in tqdm(enumerate(bedFile), total=totalEntries):
            cov = getCoverage(bamfile, entry)
            coverage[i, num] = cov
            if args.quiet:
                if i == totalEntries-1:
                    print("\rProcessing BED entry %d of %d\n100%% complete!" % (i+1, totalEntries), file=args.err)
                    args.err.flush()
                else:
                    print("\rProcessing BED entry %d of %d" %(i+1, totalEntries), end='', file=args.err)
#                               if i % checkPoint == 0:
#                                       print("\r%.1f%% complete" % (i/(checkPoint*10)), end='', file=args.err)
#                                       args.err.flush()
#                               elif i == totalEntries-1:
#                                       print("\r100% complete!", file=args.err)
#                                       args.err.flush()

args = getArgs()
mappedReads = {}
fileNames = []
totalEntries = rawLineCount(args.bed)
coverage = np.zeros((totalEntries, len(args.bams)))

if args.verbose:
    print("\n---Script started at %s---" % (str(datetime.now())), file=args.err)
    print("Total number of BED entries: %d\nTotal number of BAM files: %d\n" % (totalEntries, len(args.bams)), file=args.err)

if not args.quiet:
    print("Retrieving coverage from BAM files..", file=args.err)

for num, file in tqdm(enumerate(args.bams), total=len(args.bams), leave=True):
    processBAMFile(file, num)
if not args.quiet:
    print("Done!", file=args.err)

if not args.no_header:
    printHeader(args.out, args.bed)

# If user wants raw coverage without normalization
if args.no_norm:
    with open(args.bed, 'r') as bedFile:
        for i, line in enumerate(bedFile):
            l = line.strip().split('\t')
            print(line.strip(), *coverage[i], sep="\t", file=args.out)

# Print normalized coverage
else:
    if args.scale is None:
        scaleFactor = getScale(mappedReads)
    else:
        scaleFactor = args.scale

    if not args.quiet:
        print("Normalizing coverage for scale of %d million" % (scaleFactor/1000000), file=args.err)

    scale = [mappedReads[file]/scaleFactor for file in fileNames]
    covArray = normalizeCoverage(coverage, scale)
    if not args.quiet:
        print("Done!", file=args.err)

    with open(args.bed, 'r') as bedFile:
        for i, line in enumerate(bedFile):
            print(line.strip(), *np.round(covArray[i], 2), sep="\t", file=args.out)

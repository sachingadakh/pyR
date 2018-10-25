#! /usr/bin/python3
import sys
from random import randint

usage = 'Usage: python random_in_range.py <start number> <end number> <number of elements needed>'

if len(sys.argv) < 4:
    print(usage)
    sys.exit()

start = int(sys.argv[1])
end = int(sys.argv[2])
num = int(sys.argv[3])

for i in range(num):
    print(randint(start, end), file=sys.stdout)

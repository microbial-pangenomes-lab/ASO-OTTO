#!/bin/env python

import sys

taxdict={}
with open(sys.argv[1], "r") as taxaranks:
	header=taxaranks.readline()
	for line in taxaranks:
		splili=line.rstrip().split("\t")
		taxdict[f"{splili[0]}"]='\t'.join(splili[1::])

with open(sys.argv[3], "w") as outfile:
	outfile.write(f"sample_ID\t{header}")
	with open(sys.argv[2], "r") as samp_spec:
		for line in samp_spec:
			splili=line.rstrip().split(" ")
			tax=" ".join(splili[1::])
			outfile.write(f"{splili[0]}\t{tax}\t{taxdict[tax]}\n")

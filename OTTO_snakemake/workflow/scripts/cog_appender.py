#!/bin/env python

import sys
from io import StringIO as sio

fasta = sys.argv[1]
egg_annot = sys.argv[2]

COGdict = {}
with open(f"{egg_annot}","r") as COGfi:
	for line in COGfi:
		if line.startswith("#"): #ignore comments
			continue
		#use "gene::sample" as key and COG as item in COG-dictionary
		COGdict[line.split("\t")[0].lstrip().rstrip()]=line.split("\t")[4].split("@")[0]

buff=sio()
with open(f"{fasta}", "r") as fa:
	for line in fa:
		if line.startswith(">"):
			#use fasta header (>gene::sample) without ">" as key for COGdict to get COG and add it in the header
			buff.write(f"{line.rstrip()}::{COGdict.get(line.lstrip('>').rstrip(),'')}\n")
		else:
			buff.write(line)

sys.stdout.write(buff.getvalue())

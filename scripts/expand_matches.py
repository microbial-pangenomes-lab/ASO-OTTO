#!/usr/bin/env python

import os
import sys
import gzip
import hashlib
import subprocess
from multiprocessing import Pool

# takes centroid indeces (15 per core) and greps them 
# from files in the hash tree-like folder system
def process(item):
    d = {}
    mfile, hashes = item
    shashes = sorted(hashes)
    for i in range(0, len(shashes), 15):
        bhashes = ' -e '.join(shashes[i: i + 15])
        zgrep = subprocess.Popen(f'zgrep -e {bhashes} {mfile}',
                                 stdout=subprocess.PIPE,
                                 shell=True, universal_newlines=True)
        for k in zgrep.stdout:
            chash, member = k.rstrip().split('\t')
            d[chash] = d.get(chash, set())
            d[chash].add(member)
    return d

# helper function to paralellize the result expansion
def parallelize(inputs, cores=1):
    with Pool(processes=cores) as pool:
        return pool.map(process, inputs)

def main(indir, cores, out):
    # split input into ASO sequence and centroid indeces
    s2m={}
    f2m={}
    for l in sys.stdin:
        s = l.rstrip().split(',')
        if len(s) < 2:
            continue
        sequence = s[0]
        s2m[sequence] = s2m.get(sequence, set())

        for cha5 in s[1:]:
            if cha5 != "":
                #add hit centroid IDs to set
                s2m[sequence].add(cha5)
            # create path to desired file by slicing the five character index
            fname = f'{indir}/{cha5[:2]}/{cha5[:3]}/{cha5[:4]}.txt.gz'
            f2m[fname] = f2m.get(fname, set())
            f2m[fname].add(cha5)
    result_li = parallelize(f2m.items(), cores=cores)
    result_di = {}
    # convert results to dict
    for setdict in result_li:
        result_di.update(setdict)
    sys.stdout.write("sample_ID,gene_ID,ASO,mismatch,orig_idx,COG\n")
    # write to file
    for sequence, cha5s in s2m.items():
        spliseq = sequence.split(":")
        for cha5 in cha5s:
            if cha5 in result_di:
                for gene in result_di[cha5]:
                    splige = gene.split("::")
                    sys.stdout.write(f"{splige[1]},{splige[0]},{spliseq[1]},{spliseq[0]},{spliseq[2]},{splige[2]}\n")

if __name__ == "__main__":
    indir = sys.argv[1]
    cores = int(sys.argv[2])
    outdir = sys.argv[3]
    main(indir, cores, outdir)

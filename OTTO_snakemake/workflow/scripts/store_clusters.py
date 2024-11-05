#!/bin/env python


import os
import sys
import gzip
import hashlib


def parse_fasta(input_stream):
    s = ''
    sid = ''
    for l in input_stream:
        if l.startswith('>') and s != '':
            yield sid, s
            s = ''
            sid = l.rstrip()[1:]
        elif l.startswith('>'):
            sid = l.rstrip()[1:]
            continue
        else:
            s += l.rstrip()
    if s != '':
        yield sid, s


def main():
    clusters_file = sys.argv[1]
    out_dir = sys.argv[2]
########################################################INTERESTING PART#####################################################
    seqs = {}
    with open(clusters_file, 'r') as f: # Iterating through cluster-fasta. we can give these arbitrary headers (e.g. numbers from 00000001-60000000)
        for i, (sid, seq) in enumerate(parse_fasta(f)):
            if not i % 1_000_000:
                i1 = '%d' % (i / 1_000_000)
                sys.stderr.write(f'Processed: {i1}M clusters\n')
            # hashed = hashlib.sha256(sid.encode()).hexdigest() # can theoretically be skipped
            seqs[seq] = sid # dictionary with sequence to clusterID relationship (does not matter if the ID is hashed or not)
        i1 = '%.2f' % (i / 1_000_000)
        sys.stderr.write(f'Finished: {i1}M clusters\n')
    # centroids = set(seqs.values()) # unused line?

    tmp = {}
    for i, (sid, seq) in enumerate(parse_fasta(sys.stdin)):
        if not i % 1_000_000:
            i1 = '%d' % (i / 1_000_000)
            sys.stderr.write(f'Processed: {i1}M sequences\n')
        if seq not in seqs:
            continue
        chash = seqs[seq] # gets the hashed clusterID for the fitting sequence (means that hashed clusterIDs are correlated to sequenceIDs based on sequence)
        if chash not in tmp:
            tmp[chash] = set()
        tmp[chash].add(sid) # stores all sequenceIDs corresponding to a hashed clusterID

        if i % 100_000_000:
            continue

        if i == 0:
            continue
#########################################we'd have to rewrite this part to accomodate the new clusterIDs (if we want to go with numbers from 00000001-60000000)
        sys.stderr.write('Checkpoint: writing to files\n')
        for chash, members in tmp.items():
            folder1 = chash[:2]
            folder2 = chash[:3]
            try:
                os.mkdir(os.path.join(out_dir, folder1))
            except FileExistsError:
                pass
            try:
                os.mkdir(os.path.join(out_dir, folder1, folder2))
            except FileExistsError:
                pass
            short_hash = chash[:4]
            ofile = gzip.open(os.path.join(out_dir, folder1, folder2, f'{short_hash}.txt.gz'),
                              'at')
            for sid in members:
                ofile.write(f'{chash}\t{sid}\n')
            ofile.close()
        
        tmp = {}
    if tmp:
        sys.stderr.write('Final write: writing remaining sequences to files\n')
        for chash, members in tmp.items():
            folder1 = chash[:2]
            folder2 = chash[:3]
            try:
                os.mkdir(os.path.join(out_dir, folder1))
            except FileExistsError:
                pass
            try:
                os.mkdir(os.path.join(out_dir, folder1, folder2))
            except FileExistsError:
                pass
            short_hash = chash[:4]
            ofile = gzip.open(os.path.join(out_dir, folder1, folder2, f'{short_hash}.txt.gz'), 'at')
            for sid in members:
                ofile.write(f'{chash}\t{sid}\n')
            ofile.close()
    i1 = '%.2f' % (i / 1_000_000)
    sys.stderr.write(f'Finished: {i1}M sequences\n')


if __name__ == '__main__':
    main()

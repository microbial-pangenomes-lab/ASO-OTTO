#!/bin/env python

import sys
import itertools

# iterate through input fasta and yield a (sequence-ID, sequence) pair
def parse_fasta():
    s = ''
    sid = ''
    for l in sys.stdin:
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

#generate base62 IDs (916,132,832 in total) to be given to centroids. should be enough for most datasets :)
def generate_all_ids():
    base62_chars = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
    all_ids = (''.join(comb) for comb in itertools.product(base62_chars, repeat=5))
    return all_ids

def main():
    seqs = set()
    # create generator item that holds all base62 IDs
    id_generator = generate_all_ids()
    for i, (sid, seq) in enumerate(parse_fasta()):
        # new cluster
        if seq not in seqs:
            seqs.add(seq)
            print(f'>{next(id_generator)}\n{seq}')
	#periodically give update on how many centroid sequences were gathered
        if not i % 1_000_000:
            i1 = '%d' % (i / 1_000_000)
            i2 = '%.2f' % (len(seqs) / 1_000_000)
            sys.stderr.write(f'Processed: {i1}M ({i2}M)\n')

    i1 = '%.2f' % (i / 1_000_000)
    i2 = '%.2f' % (len(seqs) / 1_000_000)
    sys.stderr.write(f'Finished: {i1}M ({i2}M)\n')


if __name__ == '__main__':
    main()

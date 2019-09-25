#!/usr/bin/env python

from __future__ import print_function
import sys
from gzopen import gzopen

try:
    range = xrange
except NameError:
    pass

def parse_fq(filename):
    """
    fastq parser
    """

    state = 0
    label = None
    qual = None

    with gzopen(filename) as f:
        for line in f:
            line = line.rstrip("\r\n")
            if len(line) == 0:
                continue

            # fastq
            if line[0] == "@":
                if label is not None:
                    yield label, "".join(seq), "".join(qual)
                state = 1
                label = line
                seq = []
                qual = []
                continue
            elif line[0] == "+":
                state = 2
                continue

            if state == 1:
                seq.append(line)
            elif state == 2:
                qual.append(line)

    if label is not None:
        yield label, "".join(seq), "".join(qual)

def kmers(seq, k):
    """
    yield all kmers in seq
    """

    assert(k > 0)
    assert(k <= len(seq))

    for i in range(len(seq)-k+1):
        kmer = seq[i:i+k]
        yield kmer

def deletions(seq):
    """
    yield all sequences with 1bp deletion
    """
    for i in range(len(seq)):
        seq_del = seq[:i] + seq[i+1:]
        yield seq_del

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("usage: {} r1.fq r2.fq".format(sys.argv[0]))
        exit(1)

    fn1 = sys.argv[1]
    fn2 = sys.argv[2]

    hairpin = "ACGCCGGCGGCAAGTGAAGCCGCCGGCGT"
    rhairpin = "ACGCCGGCGGCTTCACTTGCCGCCGGCGT"
    hps = {hp:i for i, hp in enumerate(deletions(hairpin))}
    rhps = {rhp:len(rhairpin)-i-1 for i, rhp in enumerate(deletions(rhairpin))}

    #print(hps)
    #print(rhps)

    fp1 = parse_fq(fn1)
    fp2 = parse_fq(fn2)

    while True:
        try:
            l1, s1, q1 = next(fp1)
            l2, s2, q2 = next(fp2)
        except StopIteration:
            break

        i1 = -1
        for kmer in kmers(s1, len(hairpin)-1):
            if kmer in hps:
                i1 = hps[kmer]
                break

        if i1 == 0 or i1 == len(hairpin)-1:
            continue

        i2 = -1
        for kmer in kmers(s2, len(rhairpin)-1):
            if kmer in rhps:
                i2 = rhps[kmer]
                break

        if i2 == 0 or i2 == len(rhairpin)-1:
            continue

        if i1 == -1 and i2 == -1:
            continue


        if i1 == -1:
            print(i2)
        elif i2 == -1:
            print(i1)
        elif i1 == i2:
            print(i1)
        else:
            # discordant
            #print(i1, i2)
            pass

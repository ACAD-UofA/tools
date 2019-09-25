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

def kmers_1mm(seq, k):
    """
    yield all kmers in seq and all kmers with 1 mismatch
    """

    assert(k > 0)
    assert(k <= len(seq))

    for i in range(len(seq)):
        # 1bp deletion
        seq_del = seq[:i] + seq[i+1:]
        for kmer in kmers(seq_del,k):
            yield kmer

        for j in "ACGT":
            # 1bp mismatch
            seq_mm = seq[:i] + j + kmer[i+1:]
            for kmer in kmers(seq_mm,k):
                yield kmer

            # 1bp insertion, not seen in real data
            #seq_ins = seq[:i] + j + seq[i:]
            #for kmer in kmers(seq_ins,k):
            #    yield kmer

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("usage: {} k r1.fq r2.fq".format(sys.argv[0]))
        exit(1)

    try:
        k = int(sys.argv[1])
    except ValueError:
        print("k=`{}' invalid".format(sys.argv[1]))

    fn1 = sys.argv[2]
    fn2 = sys.argv[3]

    hairpin = "ACGCCGGCGGCAAGTGAAGCCGCCGGCGT"
    rhairpin = "ACGCCGGCGGCTTCACTTGCCGCCGGCGT"
    hpkm = frozenset(kmers_1mm(hairpin, k))
    rhpkm = frozenset(kmers_1mm(rhairpin, k))

    fp1 = parse_fq(fn1)
    fp2 = parse_fq(fn2)

    while True:
        try:
            l1, s1, q1 = next(fp1)
            l2, s2, q2 = next(fp2)
        except StopIteration:
            break

        n1 = n2 = 0
        for kmer in kmers(s1, k):
            if kmer in hpkm:
                n1 += 1
        for kmer in kmers(s2, k):
            if kmer in rhpkm:
                n2 += 1

        print(l1.split()[0], max(n1,n2), sep="\t")

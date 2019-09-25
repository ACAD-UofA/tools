#!/usr/bin/env python

from __future__ import print_function
import sys
from parsers import parse_fa
import collections
import numpy as np

def print_seq(contig_header, seq, linelen=60):
    print(contig_header)
    i = 0
    while i<len(seq):
        print(seq[i:i+linelen])
        i += linelen

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("usage: {} file.fa".format(sys.argv[0]), file=sys.stderr)
        exit(1)

    for label, seq in parse_fa(sys.argv[1]):
        c = collections.Counter(seq)
        unmasked = np.sum(c[x] for x in "ACGT")
        if unmasked == 121:
            print_seq(label, seq)

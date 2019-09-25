#!/usr/bin/env python

from __future__ import print_function
import sys

try:
    range = xrange
except NameError:
    pass

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("usage: {} k seq".format(sys.argv[0]))
        exit(1)

    try:
        k = int(sys.argv[1])
    except ValueError:
        print("k=`{}' invalid".format(sys.argv[1]))

    seq = sys.argv[2]

    kmers = set()

    for i in range(len(seq)-k):
        kmers.add(seq[i:i+k])

    print(len(kmers))
    print(kmers)

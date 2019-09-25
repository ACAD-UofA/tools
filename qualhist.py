#!/usr/bin/env python

from __future__ import print_function
import sys
import parsers
import collections

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(len(sys.argv), sys.argv)
        print("usage: {} in.fq".format(sys.argv[0]), file=sys.stderr)
        exit(1)

    quals = collections.Counter()

    for l, s, q in parsers.parse_fa(sys.argv[1]):
        quals.update(q)

    for k in sorted(quals.keys(), key=ord):
        v = quals[k]
        print("``{}''".format(k), ord(k), v)

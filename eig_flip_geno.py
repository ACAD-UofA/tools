#!/usr/bin/env python
"""
From Eigensoft 6.0.1 CONVERTF/README:
``The genotype file contains 1 line per SNP.  
  Each line contains 1 character per individual:
  0 means zero copies of reference allele.
  1 means one copy of reference allele.
  2 means two copies of reference allele.
  9 means missing data.''

But 0/2 appear to be switched in some data files. Flip them back.

Graham Gower, 2015
"""

import sys

def flip_geno(geno_file):
    with open(geno_file) as f:
        for line in f:
            for c in line:
                if c == "2":
                    c = "0"
                elif c == "0":
                    c = "2"
                sys.stdout.write(c)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("usage: {} input.geno".format(sys.argv[0]))
        sys.exit(1)

    flip_geno(sys.argv[1])

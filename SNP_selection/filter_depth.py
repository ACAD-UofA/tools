#!/usr/bin/env python
# parse vcfs and output depth histogram

from __future__ import print_function
import sys
from parsers import vcf_parser
import numpy as np
import collections

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("usage: {} file.vcf DP1 [...DPN]".format(sys.argv[0]), file=sys.stderr)
        exit(1)

    dp_match = set(map(int, sys.argv[2:]))

    vp = vcf_parser(sys.argv[1], yield_samples=True, parse_genotypes=True)
    samples = next(vp)
    if len(samples) != 1:
        raise Exception("{} samples. Single sample vcf expected.".format(len(samples)))

    s = samples[0]

    for vline in vp:
        chrom = vline[2]
        pos = vline[3]
        qual = float(vline[7])
        genotype_field = vline[11][s]
        ad = genotype_field["AD"]
        dp = np.sum(map(int, ad.split(",")))
        if dp in dp_match:
            print(chrom, pos-1, pos, sep="\t")

#!/usr/bin/env python
# parse vcfs and output depth histogram

from __future__ import print_function
import sys
from parsers import vcf_parser
import numpy as np
import collections

def vcf_dp(fn_list, max_dp=1024):
    dp_array = np.zeros(max_dp, dtype=int)
    qual_dictlist = collections.defaultdict(list)
    for filename in fn_list:
        vp = vcf_parser(filename, yield_samples=True, parse_genotypes=True)
        samples = next(vp)
        if len(samples) != 1:
            raise Exception("{} samples. Single sample vcf expected.".format(len(samples)))
        s = samples[0]
        for vline in vp:
            qual = float(vline[7])
            genotype_field = vline[11][s]
            ad = genotype_field["AD"]
            dp = np.sum(map(int, ad.split(",")))
            if dp < max_dp:
                dp_array[dp] += 1
                qual_dictlist[dp].append(qual)

    qual_mean = [np.mean(qual_dictlist[i]) if len(qual_dictlist[i])>0 else 0 for i in range(max_dp)]
    qual_var = [np.var(qual_dictlist[i]) if len(qual_dictlist[i])>0 else 0 for i in range(max_dp)]
    return dp_array, qual_mean, qual_var

def main():
    dp_array, qmean, qvar = vcf_dp(sys.argv[1:])
    tot = np.sum(dp_array)

    for dp, (n, q, qv) in enumerate(zip(dp_array, qmean, qvar)):
        print(dp, n, float(n)/tot, q, qv, sep="\t")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise Exception("usage: {} file1.vcf [...fileN.vcf]".format(sys.argv[0]))

    try:
        main()
    except KeyboardInterrupt:
        pass


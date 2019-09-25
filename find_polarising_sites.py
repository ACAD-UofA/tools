#!/usr/bin/env python

from __future__ import print_function
import sys
import re

def parse_geno(filename):
    with open(filename) as f:
        for line in f:
            yield line

def parse_ind(filename, a, b):

    re_a = [re.compile(x) for x in a]
    re_b = [re.compile(x) for x in b]

    ai, bi = set(), set()
    indlist = []

    with open(filename) as f:
        for i, line in enumerate(f):
            fields = line.split()
            ind = fields[0]
            indlist.append(ind)
            #sex = fields[1]
            group = fields[2]
            for r in re_a:
                if r.search(group):
                    ai.add(i)
                    break
            for r in re_b:
                if r.search(group):
                    bi.add(i)
                    break

    both = ai&bi
    if len(both) > 0:
        for i in both:
            print("Error: {} in group a and group b".format(indlist[i]), file=sys.stderr)
        exit(1)

    return ai, bi, indlist

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="count SNPs with one allele fixed in group a, different allele fixed in group b)")
    parser.add_argument("--print-groups", action="store_true", help="print individuals for each group")
    parser.add_argument("--print-sites", action="store_true", help="print SNP number for each polarising site")
    parser.add_argument("-a", "--pop-a", required=True, action="append", help="regex for group a")
    parser.add_argument("-b", "--pop-b", required=True, action="append", help="regex for group b")
    parser.add_argument("geno_fn", metavar="fn.geno", help="eigenstrat genotypes")
    parser.add_argument("ind_fn", metavar="fn.ind", help="eigenstrat individuals")
    #parser.add_argument("snp_fn", metavar="fn.snp", help="eigenstrat SNP list")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    ai, bi, indlist = parse_ind(args.ind_fn, args.pop_a, args.pop_b)

    # number of 'a' and 'b' individuals
    na = len(ai)
    nb = len(bi)

    if args.print_groups:
        print("group a\n--------")
        for i in ai:
            print(indlist[i])

        print("group b\n--------")
        for i in bi:
            print(indlist[i])

    nsites = npolar = 0

    for lineno, line in enumerate(parse_geno(args.geno_fn), 1):

        line = line.strip()

        a = b = 0
        a_obs = b_obs = 0

        for i, gt in enumerate(line):
            if gt == "9":
                continue
            if i in ai:
                a += int(gt)
                a_obs += 1
            elif i in bi:
                b += int(gt)
                b_obs += 1

        if a_obs < 0.5*na or b_obs < 0.5*nb:
            continue

        nsites += 1

        if ((a == 2*a_obs and b == 0) or
                (b == 2*b_obs and a == 0)):
            npolar += 1
            if args.print_sites:
                print(lineno)

    print(npolar, "polarising sites out of", nsites, "comparable sites")

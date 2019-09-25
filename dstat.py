#!/usr/bin/env python

from __future__ import print_function
import re
import sys
import collections
import itertools
import random
import csv
import numpy as np

def parse_lines(filename):
    with open(filename) as f:
        for line in f:
            yield line

def mk_eig_chr_map_inv():
    """
    Undo eigensoft's awkward chromosome naming convention.
    """
    chrmap = {"23":"X", "24":"Y", "90":"M"}
    for i in range(23):
        chrmap[str(i)] = str(i)
    for i in range(33, 90):
        chrmap[str(i)] = str(i-10)
    return chrmap

def parse_snp(filename):
    re_s = re.compile(r"[ \t]+")
    chrmap = mk_eig_chr_map_inv()
    lineno = 0
    snps = {}
    for line in parse_lines(filename):
        fields = re_s.split(line.strip())

        snpid = fields[0]
        chrom = chrmap[fields[1]]
        gpos = fields[2]
        pos = int(fields[3])

        snps[lineno] = (snpid, chrom, gpos, pos)
        lineno += 1

    return snps

def parse_ind(filename):
    re_s = re.compile(r"[ \t]+")
    lineno = 0
    groups = collections.defaultdict(list)
    for line in parse_lines(filename):
        fields = re_s.split(line.strip())

        #ind = fields[0]
        #sex = fields[1]
        group = fields[2]

        groups[group].append(lineno)
        lineno += 1

    return groups

def dstat(genotypes, groups):
    """
    Calculate D statistic. See Patterson et al. (2012).
    """

    skip = False
    af = {}
    numerator = 0
    denominator = 0
    lineno = 0

    baba = []
    abba = []
    covered = []

    for line in genotypes:
        lineno += 1

        # Find allele frequency for each group.
        for i, g in enumerate(groups):
            gts = [int(line[j]) for j in gind[g] if line[j]!="9"]
            if len(gts) == 0:
                # No observations for this site.
                skip = True
                break
            # Frequency of the major allele.
            af[i] = float(sum(gts)) / (2*len(gts))

        if skip:
            skip = False
            continue
        
        covered.append(lineno-1)

        w = af[0]
        x = af[1]
        y = af[2]
        z = af[3]

        # Follow the BABA-ABBA convention, as in AdmixTools.
        num = (w - x)*(y - z)
        den = (w + x - 2*w*x)*(y + z - 2*y*z)
        numerator += num
        denominator += den

        if num > 0:
            baba.append(lineno-1)
        elif num < 0:
            abba.append(lineno-1)

    try:
        d = numerator / denominator
    except ZeroDivisionError:
        d = 0

    return d, frozenset(baba), frozenset(abba), frozenset(covered)


def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="Calculate D statistic from *.{snp,geno} files.")
    parser.add_argument("--snp", required=True, help="eigensoft .snp file")
    parser.add_argument("--ind", required=True, help="eigensoft .ind file")
    parser.add_argument("--geno", required=True, help="eigensoft .geno file")
    parser.add_argument("--groups", nargs='+', required=True, help="P1,P2,P3,O")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    gind = parse_ind(args.ind)
    groupslist = [groups.split(",") for groups in args.groups]

    for groups in groupslist:
        if len(groups) != 4:
            print("Error: Must specify 4 groups", file=sys.stderr)
            exit(1)

        for g in groups:
            if g not in gind:
                print("Error: group ``{}'' not in {}.".format(g, args.ind), file=sys.stderr)
                exit(1)

    snps = parse_snp(args.snp)
    genotypes = list(parse_lines(args.geno))

    dd = {}
    for groups, groups_str in zip(groupslist, args.groups):
        dd[groups_str] = dstat(genotypes, groups)

    writer = csv.DictWriter(sys.stdout, delimiter="\t", fieldnames=[
        "TEST1", "TEST2",
        "TEST1 covered SNPs", "TEST1 informative SNPs",
        "TEST2 covered SNPs", "TEST2 informative SNPs",
        "BABA1",
        "ABBA1",
        "BABA2",
        "ABBA2",
        "BABA agree", "ABBA agree",
        "disagree", "uninformative"])
    writer.writeheader()

    for i, (g1, g1_str) in enumerate(zip(groupslist, args.groups),1):
        for g2, g2_str in zip(groupslist[i:], args.groups[i:]):

            d1, baba1, abba1, covered1 = dd[g1_str]
            d2, baba2, abba2, covered2 = dd[g2_str]

            assert(len(baba1 & abba1) == 0)
            assert(len(baba2 & abba2) == 0)

            baba_agree = baba1 & baba2
            abba_agree = abba1 & abba2
            disagree = (baba1 & abba2) | (baba2 & abba1)

            informative1 = baba1 | abba1
            informative2 = baba2 | abba2
            uninformative = (covered1 - informative1) | (covered2 - informative2)

            baba_uninformative = (baba1 - baba2 - abba2) | (baba2 - baba1 - abba1)
            abba_uninformative = (abba1 - abba2 - baba2) | (abba2 - abba1 - baba1)
            partial_uninformative = baba_uninformative | abba_uninformative

            row = {}
            row["TEST1"] = "((({}, {}), {}), {})".format(*g1)
            row["TEST2"] = "((({}, {}), {}), {})".format(*g2)

            row["TEST1 covered SNPs"] = len(covered1)
            row["TEST1 informative SNPs"] = len(informative1)
            row["TEST2 covered SNPs"] = len(covered2)
            row["TEST2 informative SNPs"] = len(informative2)

            row["disagree"] = len(disagree)
            row["uninformative"] = len(partial_uninformative)

            row["BABA1"] = len(baba1)
            row["BABA2"] = len(baba2)
            row["BABA agree"] = len(baba_agree)

            row["ABBA1"] = len(abba1)
            row["ABBA2"] = len(abba2)
            row["ABBA agree"] = len(abba_agree)

            writer.writerow(row)

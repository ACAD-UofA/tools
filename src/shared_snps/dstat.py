#!/usr/bin/env python

from __future__ import print_function
import re
import sys
import collections
import csv
from scipy.stats import hypergeom

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

def dstat(genotypes, groups, gind):
    """
    Calculate D statistic. See Patterson et al. (2012).
    """

    num = 0
    den = 0
    snpnum = -1

    baba = []
    abba = []

    def af(gt, group, gind):
        """
        Frequency of the major allele.
        """
        gts = [int(gt[i]) for i in gind[group] if gt[i] != "9"]
        n_inds = len(gts)
        if n_inds == 0:
            return None
        return float(sum(gts)) / n_inds

    for row in genotypes:
        snpnum += 1

        w = af(row, groups[0], gind)
        x = af(row, groups[1], gind)
        y = af(row, groups[2], gind)
        z = af(row, groups[3], gind)

        if None in (w, x, y, z):
            continue

        # Follow the BABA-ABBA convention, as in AdmixTools.
        num += (w - x)*(y - z)
        den += (w + x - 2*w*x)*(y + z - 2*y*z)

        if (abs(w-y) < abs(w-x) and abs(w-y) < abs(y-x) and
            abs(w-y) < abs(w-z) and abs(w-y) < abs(y-z) and
            abs(x-z) < abs(w-x) and abs(x-z) < abs(y-x) and
            abs(x-z) < abs(w-z) and abs(x-z) < abs(y-z)):
            baba.append(snpnum)
        elif (abs(x-y) < abs(x-w) and abs(x-y) < abs(y-w) and
              abs(x-y) < abs(x-z) and abs(x-y) < abs(y-z) and
              abs(w-z) < abs(x-w) and abs(w-z) < abs(y-w) and
              abs(w-z) < abs(x-z) and abs(w-z) < abs(y-z)):
            abba.append(snpnum)

    try:
        d = num / den
    except ZeroDivisionError:
        d = 0

    return d, frozenset(baba), frozenset(abba)


def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="Calculate D statistic from *.{snp,geno} files.")
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

    genotypes = list(parse_lines(args.geno))

    dd = {}
    for groups, groups_str in zip(groupslist, args.groups):
        dd[groups_str] = dstat(genotypes, groups, gind)

    writer = csv.DictWriter(sys.stdout, delimiter="\t", fieldnames=[
        "1", "2",
        "BABA | ABBA",
        "BABA1",
        "ABBA1",
        "BABA2",
        "ABBA2",
        "BABA1 & BABA2",
        "BABA1 & ABBA2",
        "ABBA1 & BABA2",
        "ABBA1 & ABBA2",
        "p"
        ])
    writer.writeheader()

    for i, (g1, g1_str) in enumerate(zip(groupslist, args.groups),1):
        for g2, g2_str in zip(groupslist[i:], args.groups[i:]):

            d1, baba1, abba1 = dd[g1_str]
            d2, baba2, abba2 = dd[g2_str]

            # keep only those sites for which we have complete data
            informative = (baba1 | abba1) & (baba2 | abba2)

            row = {}
            row["1"] = "((({}, {}), {}), {})".format(*g1)
            row["2"] = "((({}, {}), {}), {})".format(*g2)

            row["BABA | ABBA"] = len(informative)
            baba1 &= informative
            abba1 &= informative
            baba2 &= informative
            abba2 &= informative

            row["BABA1"] = len(baba1)
            row["BABA2"] = len(baba2)
            row["ABBA1"] = len(abba1)
            row["ABBA2"] = len(abba2)

            row["BABA1 & BABA2"] = len(baba1 & baba2)
            row["BABA1 & ABBA2"] = len(baba1 & abba2)
            row["ABBA1 & BABA2"] = len(abba1 & baba2)
            row["ABBA1 & ABBA2"] = len(abba1 & abba2)

            # chance of observing len(baba1&baba2) or more sites in common
            #p = 1 -hypergeom.cdf(len(baba1&baba2)-1, len(informative), len(baba1), len(baba2))
            p = hypergeom.sf(len(baba1&baba2)-1, len(informative), len(baba1), len(baba2))
            row["p"] = p

            writer.writerow(row)

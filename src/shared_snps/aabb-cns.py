#!/usr/bin/env python

from __future__ import print_function
import re
import sys
import collections
import csv
from scipy.stats import hypergeom
from scipy import random

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

def aabb_count_old(genotypes, groups, gind):
    """
    Identify shared derived SNPs.
    
    genotypes -- rows of genotypes from the .geno file
    groups    -- list of 4 groups, (W, X, Y, O)
    gind      -- map from the group name to a list of genotype columns

    The group W is assumed to have a hybrid ancestry, with unknown
    proportions of both X and Y.  O is a more distantly related outgroup.
    We wish to identify derived alleles from X which are also found in W.
    The set of candidate sites are observed to have the pattern (?,A,B,B),
    where the shared derived sites have the pattern (A,A,B,B).
    """

    snpnum = -1
    aabb = []
    babb = []

    def af(gt, group, gind):
        """Calculate allele frequency."""
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

        # mandate (?,A,B,B)
        if abs(z-y) >= abs(z-x) or abs(z-y) >= abs(y-x):
            continue

        if (abs(z-y) < abs(z-w) and abs(z-y) < abs(y-w) and
                abs(w-x) < abs(w-y) and abs(w-x) < abs(w-z) and
                abs(w-x) < abs(x-y) and abs(w-x) < abs(x-y)):
            # (A,A,B,B)
            aabb.append(snpnum)
        else:
            # (B,A,B,B)
            babb.append(snpnum)

    return frozenset(aabb), frozenset(babb)

def aabb_count(genotypes, groups, gind):
    """
    Identify shared derived SNPs.
    
    genotypes -- rows of genotypes from the .geno file
    groups    -- list of 4 groups, (W, X, Y, O)
    gind      -- map from the group name to a list of genotype columns

    The group W is assumed to have a hybrid ancestry, with unknown
    proportions of both X and Y.  O is a more distantly related outgroup.
    We wish to identify derived alleles from X which are also found in W.
    The set of candidate sites are observed to have the pattern (?,A,B,B),
    where the shared derived sites have the pattern (A,A,B,B).
    """

    snpnum = -1
    aabb = []
    babb = []

    def cns(gt, group, gind):
        """Return consensus allele, or randomly choose if no consensus."""
        gts = [int(gt[i]) for i in gind[group] if gt[i] != "9"]
        n_inds = len(gts)
        if n_inds == 0:
            return None

        s_gts = sum(gts)
        if s_gts > len(gts):
            # greater than 50% of alleles are REF
            return 1
        elif s_gts < len(gts):
            # fewer than 50% of alleles are REF
            return 0
        else:
            return None #random.choice((0,1))

    for row in genotypes:
        snpnum += 1

        w = cns(row, groups[0], gind)
        x = cns(row, groups[1], gind)
        y = cns(row, groups[2], gind)
        z = cns(row, groups[3], gind)

        if None in (w, x, y, z):
            continue

        # mandate (?,A,B,B)
        if y!=z or x == y:
            continue

        if w == x:
            # (A,A,B,B)
            aabb.append(snpnum)
        else:
            # (B,A,B,B)
            babb.append(snpnum)

    return frozenset(aabb), frozenset(babb)

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="Identify shared derived SNPs from *.{snp,geno} files.")
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

    ab = {}
    for groups, groups_str in zip(groupslist, args.groups):
        ab[groups_str] = aabb_count(genotypes, groups, gind)

    writer = csv.DictWriter(sys.stdout, delimiter="\t", fieldnames=[
        "1", "2",
        "?ABB",
        "AABB1", "BABB1",
        "AABB2", "BABB2",
        "AABB1 & AABB2",
        "AABB1 & BABB2", "AABB2 & BABB1",
        "BABB1 & BABB2",
        "p"
        ])
        
    writer.writeheader()

    for i, (g1, g1_str) in enumerate(zip(groupslist, args.groups),1):
        for g2, g2_str in zip(groupslist[i:], args.groups[i:]):

            aabb1, babb1 = ab[g1_str]
            aabb2, babb2 = ab[g2_str]

            # keep only those sites for which we have complete data
            informative = (aabb1|babb1) & (aabb2|babb2)
            aabb1 &= informative
            aabb2 &= informative
            babb1 &= informative
            babb2 &= informative

            row = {}
            row["1"] = "({}, {}, {}, {})".format(*g1)
            row["2"] = "({}, {}, {}, {})".format(*g2)

            row["?ABB"] = len(informative)

            row["AABB1"] = len(aabb1)
            row["AABB2"] = len(aabb2)
            row["BABB1"] = len(babb1)
            row["BABB2"] = len(babb2)

            row["AABB1 & AABB2"] = len(aabb1&aabb2)
            row["AABB1 & BABB2"] = len(aabb1&babb2)
            row["AABB2 & BABB1"] = len(aabb2&babb1)
            row["BABB1 & BABB2"] = len(babb1&babb2)

            # chance of observing len(aabb1&aaab2) or more sites in common
            #p = 1 -hypergeom.cdf(len(aabb1&aabb2)-1, len(informative), len(aabb1), len(aabb2))
            p = hypergeom.sf(len(aabb1&aabb2)-1, len(informative), len(aabb1), len(aabb2))
            row["p"] = p

            writer.writerow(row)

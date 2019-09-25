#!/usr/bin/env python
# Copyright (c) 2015 Graham Gower <graham.gower@gmail.com>
#
# Permission to use, copy, modify, and distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.
#
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
# ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
# ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
# OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

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

def parse_snp(filename):
    re_s = re.compile(r"[ \t]+")
    chrmap = mk_eig_chr_map_inv()
    snps = {}
    for lineno, line in enumerate(parse_lines(filename)):
        fields = re_s.split(line.strip())

        snpid = fields[0]
        chrom = chrmap[fields[1]]
        gpos = fields[2]
        pos = int(fields[3])

        snps[lineno] = (snpid, chrom, gpos, pos)

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

def aabb_count_af(genotypes, groups, gind):
    """
    Identify shared derived SNPs.
    
    genotypes -- rows of genotypes from the .geno file
    groups    -- list of 4 groups, (W, X, Y, Z)
    gind      -- map from the group name to a list of genotype columns

    The group W is assumed to have a hybrid ancestry, with unknown
    proportions of both X and Y.  Z is a more distantly related outgroup.
    We wish to identify derived alleles from X which are also found in W.
    The set of candidate sites are observed to have the pattern (?,A,B,B),
    where the shared derived sites have the pattern (A,A,B,B).
    """

    aabb = []
    babb = []

    def af(gt, group, gind):
        """Calculate allele frequency."""
        gts = [int(gt[i]) for i in gind[group] if gt[i] != "9"]
        n_inds = len(gts)
        if n_inds == 0:
            return None
        freq = float(sum(gts)) / n_inds

        n_alleles = 2*n_inds
        s_gts = sum(gts)

        # estimator of heterozygosity per allele
        het = float(s_gts * (n_alleles-s_gts)) / (n_alleles*n_alleles*(n_alleles-1))

        return freq, het

    class NegativeF2(Exception):
        pass

    def f2(a, b):
        """Unbiased F2 estimator, see Patterson et al. (2012), Appendix A."""
        f2_ = (a[0]-b[0])**2 -a[1] -b[1]
        if f2_ < 0:
            raise NegativeF2
        return f2_

    for snpnum, row in enumerate(genotypes):
        w = af(row, groups[0], gind)
        x = af(row, groups[1], gind)
        y = af(row, groups[2], gind)
        z = af(row, groups[3], gind)

        if None in (w, x, y, z):
            continue

        try:
            # mandate (?,A,B,B)
            #if abs(z-y) >= abs(z-x) or abs(z-y) >= abs(y-x):
            #    continue
            if f2(z,y) >= f2(z,x) or f2(z,y) >= f2(y,x):
                continue

            #if (abs(z-y) < abs(z-w) and abs(z-y) < abs(y-w) and
            #        abs(w-x) < abs(w-y) and abs(w-x) < abs(w-z) and
            #        abs(w-x) < abs(x-y) and abs(w-x) < abs(x-y)):
            #    # (A,A,B,B)
            #    aabb.append(snpnum)
            #else:
            #    # (B,A,B,B)
            #    babb.append(snpnum)
            if (f2(z,y) < f2(z,w) and f2(z,y) < f2(y,w) and
                    f2(w,x) < f2(w,y) and f2(w,x) < f2(w,z) and
                    f2(w,x) < f2(x,y) and f2(w,x) < f2(x,z)):
                aabb.append(snpnum)
            elif (f2(w,x) > f2(w,y) and f2(w,x) > f2(w,z) and
                    f2(x,y) > f2(w,y) and f2(x,y) > f2(w,z) and
                    f2(x,z) > f2(w,y) and f2(x,z) > f2(w,z)):
                babb.append(snpnum)

        except NegativeF2:
            continue

    return frozenset(aabb), frozenset(babb)

def aabb_count_cns(genotypes, groups, gind):
    """
    Identify shared derived SNPs.
    
    genotypes -- rows of genotypes from the .geno file
    groups    -- list of 4 groups, (W, X, Y, Z)
    gind      -- map from the group name to a list of genotype columns

    The group W is assumed to have a hybrid ancestry, with unknown
    proportions of both X and Y.  Z is a more distantly related outgroup.
    We wish to identify derived alleles from X which are also found in W.
    The set of candidate sites are observed to have the pattern (?,A,B,B),
    where the shared derived sites have the pattern (A,A,B,B).
    """

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

    for snpnum, row in enumerate(genotypes):
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
    parser.add_argument("-c", "--consensus", action="store_true", default=False, help="Use consensus allele, not allele frequency")
    parser.add_argument("-i", "--print-ids", action="store_true", default=False, help="Output common SNP ids for all comparisons separate files")
    parser.add_argument("--snp", required=True, help="eigensoft .snp file")
    parser.add_argument("--ind", required=True, help="eigensoft .ind file")
    parser.add_argument("--geno", required=True, help="eigensoft .geno file")
    parser.add_argument("--groups", nargs='+', required=True, help="P1,P2,P3,O")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    snpdict = parse_snp(args.snp)
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

    if args.consensus:
        aabb_count = aabb_count_cns
    else:
        aabb_count = aabb_count_af

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

            assert(g1[1:] == g2[1:])

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

            if not args.print_ids:
                continue

            # output the SNP ids found to be in common
            fname = "{}-SNPids_{}-vs-{}.txt".format(g1[1], g1[0], g2[0])
            with open(fname, "w") as fout:
                snpids = [snpdict[idx][0] for idx in (aabb1&aabb2)]
                snpids.sort()
                fout.write("\n".join(snpids))
                fout.write("\n")

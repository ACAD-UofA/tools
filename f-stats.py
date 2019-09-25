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

def pop_params(gt, group, gind):
    """
    Calculate population parameters from the given getotypes.
    
    Returns:
        consensus allele,
        allele frequency,
        heterozygosity,
        the number of alleles observed.
    """

    gts = [int(gt[i]) for i in gind[group] if gt[i] != "9"]
    n_inds = len(gts)
    if n_inds == 0:
        return None, None, None, None

    s_gts = sum(gts)
    if s_gts > len(gts):
        # greater than 50% of alleles are REF
        cns_allele = 1
    elif s_gts < len(gts):
        # fewer than 50% of alleles are REF
        cns_allele = 0
    else:
        cns_allele = None

    af = float(s_gts) / n_inds

    n_alleles = 2*len(gts)
    het = float(s_gts * (n_alleles-s_gts)) / (n_alleles*(n_alleles-1))

    return cns_allele, af, het, n_alleles

def F2(genotypes, gind, p1, p2):
    """
    Calculates overlap between genetic drift paths p1<->p2 and p1<->p2.

    Uses the unbiased estimator given in Patterson et al. (2012).
    """

    f2_rows = 0
    f2_sum = 0

    for row in genotypes:
        _, af_p1, het_p1, n_p1 = pop_params(row, p1, gind)
        _, af_p2, het_p2, n_p2 = pop_params(row, p2, gind)

        if None in (af_p1, af_p2):
            continue

        f2_rows += 1
        f2_sum += (af_p1 - af_p2)**2 - het_p1/n_p1 - het_p2/n_p2

    f2 = f2_sum / f2_rows

    return f2


def F3(genotypes, gind, px, p1, p2):
    """
    Calculates overlap between genetic drift paths px<->p1 and px<->p2.

    Uses the unbiased estimator given in Patterson et al. (2012).
    """

    f3_rows = 0
    f3_sum = 0

    for row in genotypes:
        _, af_px, het_px, n_px = pop_params(row, px, gind)
        _, af_p1, _, _ = pop_params(row, p1, gind)
        _, af_p2, _, _ = pop_params(row, p2, gind)

        if None in (af_p1, af_p2, af_px):
            continue

        f3_rows += 1
        f3_sum += (af_px - af_p1)*(af_px - af_p2) - het_px/n_px

    f3 = f3_sum / f3_rows

    return f3

def F4(genotypes, gind, p1, p2, p3, p4):
    """
    Calculates overlap between genetic drift paths p1<->p2 and p3<->p4.

    This is an unbiased estimator according to Durand et al. (2011).
    """

    f4_rows = 0
    f4_sum = 0
    cns4_rows = 0
    cns4_sum = 0

    for row in genotypes:
        _, af_p1, _, _ = pop_params(row, p1, gind)
        _, af_p2, _, _ = pop_params(row, p2, gind)
        _, af_p3, _, _ = pop_params(row, p3, gind)
        _, af_p4, _, _ = pop_params(row, p4, gind)

        if None in (af_p1, af_p2, af_p3, af_p4):
            continue

        f4_rows += 1
        f4_sum += (af_p1 - af_p2)*(af_p3 - af_p4)

    f4 = f4_sum / f4_rows

    return f4

def cns(genotypes, gind, p1, p2, p3, p4):

    n_rows = 0
    babb = 0
    abbb = 0
    aaba = 0
    aabb = 0

    for row in genotypes:
        cns_p1, _, _, _ = pop_params(row, p1, gind)
        cns_p2, _, _, _ = pop_params(row, p2, gind)
        cns_p3, _, _, _ = pop_params(row, p3, gind)
        cns_p4, _, _, _ = pop_params(row, p4, gind)

        if None in (cns_p1, cns_p2, cns_p3, cns_p4):
            continue

        n_rows += 1

        ab = (cns_p1, cns_p2, cns_p3, cns_p4)

        if ab == (0,1,0,0) or ab == (1,0,1,1):
            babb += 1
        elif ab == (1,0,0,0) or ab == (0,1,1,1):
            abbb += 1
        elif ab == (0,0,1,0) or ab == (1,1,0,1):
            aaba += 1
        elif ab == (0,0,1,1) or ab == (1,1,0,0):
            aabb += 1

    return babb, abbb, aaba, aabb, n_rows



def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="Identify shared derived SNPs from *.{snp,geno} files.")
    parser.add_argument("--ind", required=True, help="eigensoft .ind file")
    parser.add_argument("--geno", required=True, help="eigensoft .geno file")
    parser.add_argument("--groups", nargs='+', required=True, help="P1,P2,P3,P4")
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

    writer = csv.DictWriter(sys.stdout, delimiter="\t", fieldnames=[
        "P1", "P2", "P3", "P4",
        "F2(P1,P2)",
        "F3(P3;P1,P4)",
        "F3(P3;P2,P4)",
        "F4(P1,P3;P2,P4)",
        "BABB",
        "ABBB",
        "AABA",
        "AABB",
        "n_sites"
        ])
        
    writer.writeheader()

    for groups in groupslist:
        p1 = groups[0]
        p2 = groups[1]
        p3 = groups[2]
        p4 = groups[3]

        f2 = F2(genotypes, gind, p1, p2)
        f3_p1 = F3(genotypes, gind, p3, p1, p4)
        f3_p2 = F3(genotypes, gind, p3, p2, p4)
        f4 = F4(genotypes, gind, p1, p3, p2, p4)

        babb, abbb, aaba, aabb, n_sites = cns(genotypes, gind, p1, p3, p2, p4)

        row = {}
        row["P1"] = p1
        row["P2"] = p2
        row["P3"] = p3
        row["P4"] = p4

        row["F2(P1,P2)"] = f2
        row["F3(P3;P1,P4)"] = f3_p1
        row["F3(P3;P2,P4)"] = f3_p2
        row["F4(P1,P3;P2,P4)"] = f4
        row["BABB"] = babb
        row["ABBB"] = abbb
        row["AABA"] = aaba
        row["AABB"] = aabb
        row["n_sites"] = n_sites

        writer.writerow(row)

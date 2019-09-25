#!/usr/bin/env python
"""
Produce .tsv, containing GT stats aggregated by DP.

Author: Graham Gower, 2015
"""

from __future__ import print_function
import gzip
import collections
import sys
import os.path
import csv

class ParseError(Exception):
    pass

def parse_vcf(vcf_in, max_dp, min_qual):
    """
    Return gt_count, aggregated by DP.
    """

    samples = None
    lineno = 0

    gt_count = collections.defaultdict(collections.Counter)

    if vcf_in.endswith(".gz"):
        vcf_open = gzip.open
    else:
        vcf_open = open

    with vcf_open(vcf_in) as f:

        # parse header
        for line in f:
            lineno += 1
            if line.startswith("#CHROM"):
                # column header line
                # #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample0 [sample1 ...]
                line = line.rstrip("\r\n")
                fields = line.split("\t")
                if len(fields) < 9:
                    raise ParseError("{}: line {}: expected at least 9 columns for vcf column header.".format(vcf_in, lineno))
                samples = fields[9:]
                break

        if samples is None:
            raise ParseError("{}: no vcf column header found.".format(vcf_in))

        # parse variants
        for line in f:
            lineno += 1
            line = line.rstrip("\r\n")
            fields = line.split("\t")

            chrom = fields[0]
            pos = int(fields[1])
            id = fields[2]
            ref = fields[3]
            alt = fields[4]
            qual = float(fields[5])
            filter = fields[6]
            info_str = fields[7]
            format_str = fields[8]
            genotypes = fields[9:]

            if qual < min_qual:
                continue

            fmt_fields = format_str.split(":")

            if "GT" not in fmt_fields:
                raise ParseError("{}: line {}: no FORMAT/GT field. Genotype has not been called.".format(vcf_in, lineno))

            if "DP" not in fmt_fields:
                raise ParseError("{}: line {}: missing FORMAT/DP field. Depth is required".format(vcf_in, lineno))

            for sample, genotype_str in zip(samples, genotypes):
                gen_fields = genotype_str.split(":")
                genotype = dict()
                for fmt, gen in zip(fmt_fields, gen_fields):
                    genotype[fmt] = gen
                
                dp = int(genotype["DP"])
                if dp > max_dp:
                    dp = max_dp

                gt_str = genotype["GT"]
                if gt_str == "./.":
                    gt = 0
                else:
                    gt = sum([int(x) for x in gt_str.split("/")])

                gt_count[dp][gt] += 1

    return gt_count


def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="Print GT stats for one or more VCFs")
#    parser.add_argument("--min-depth", type=int, default=5, help="Filter calls based upon the DP field [%(default)s]")
    parser.add_argument("--min-qual", type=int, default=50, help="Filter calls based upon the QUAL field [%(default)s]")
    parser.add_argument("vcf", nargs="+", help="One or more VCFs containing genotype data (optionally gzipped)")
    return parser.parse_args()

if __name__ == "__main__":

    args = parse_args()

    header = ["Sample"]
    for x in range(1,6):
        header.extend(["Total (DP>={})".format(x), "REF/REF (DP>={})".format(x), "REF/ALT (DP>={})".format(x), "ALT/ALT (DP>={})".format(x)])

    writer = csv.DictWriter(sys.stdout, fieldnames=header, delimiter="\t")
    writer.writeheader()

    for vcf_in in args.vcf:

        base = os.path.basename(vcf_in)
        if base.endswith(".vcf.gz"):
            base = base[:-7]
        elif base.endswith(".vcf"):
            base = base[:-4]
        if base.endswith("_10k_snps"):
            base = base[:-9]

        gt_count = parse_vcf(vcf_in, max_dp=5, min_qual=args.min_qual)

        row = {"Sample": base}

        for x in range(5,0,-1):
            hom_ref = sum([gt_count[y][0] for y in range(x,10)])
            het = sum([gt_count[y][1] for y in range(x,10)])
            hom_alt = sum([gt_count[y][2] for y in range(x,10)])

            row["Total (DP>={})".format(x)] = hom_ref + het + hom_alt
            row["REF/REF (DP>={})".format(x)] = hom_ref
            row["REF/ALT (DP>={})".format(x)] = het
            row["ALT/ALT (DP>={})".format(x)] = hom_alt

        writer.writerow(row)

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
"""
Print stats on .vcf GT numbers.
"""

from __future__ import print_function
import sys
import os.path
import collections
import gzip
import csv

class ParseError(Exception):
    pass

def parse_vcf(vcf_in, min_depth, min_qual):

    samples = None
    lineno = 0

    if vcf_in.endswith(".gz"):
        vcf_open = gzip.open
    else:
        vcf_open = open
    
    gt_counts = collections.Counter({"0/0":0, "0/1":0, "1/1":0})

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

            for sample, genotype_str in zip(samples, genotypes):
                gen_fields = genotype_str.split(":")
                genotype = dict()
                for fmt, gen in zip(fmt_fields, gen_fields):
                    genotype[fmt] = gen
                
                if int(genotype["DP"]) < min_depth:
                    continue

                gt_counts[genotype["GT"]] += 1

    return gt_counts


def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="Print GT stats for one or more VCFs")
    parser.add_argument("--min-depth", type=int, default=5, help="Filter calls based upon the DP field [%(default)s]")
    parser.add_argument("--min-qual", type=int, default=50, help="Filter calls based upon the QUAL field [%(default)s]")
    parser.add_argument("vcf", nargs="+", help="One or more VCFs containing genotype data (optionally gzipped)")
    return parser.parse_args()

if __name__ == "__main__":

    args = parse_args()

    writer = csv.DictWriter(sys.stdout, fieldnames=["Sample", "Total SNPs", "0/0", "0/1", "1/1"])
    writer.writeheader()

    for filename in args.vcf:
        base = os.path.basename(filename)[:-len(".vcf.gz")]
        gt_counts = parse_vcf(filename, args.min_depth, args.min_qual)
        gt_counts["Sample"] = base
        gt_counts["Total SNPs"] = gt_counts["0/0"] + gt_counts["0/1"] + gt_counts["1/1"]
        writer.writerow(gt_counts)


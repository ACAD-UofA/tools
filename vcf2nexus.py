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
Convert one or more vcf to nexus format.
"""

from __future__ import print_function
import sys
import collections
import parsers
import random

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="Convert the SNPs from a vcf to nexus format (with 0=REF/REF, 1=REF/ALT, 2=ALT/ALT, ?=unknown).")
    parser.add_argument("-r", "--reference", action="store_true", default=False, help="Include a reference individual.")
    parser.add_argument("-n", "--nucleotide", action="store_true", default=False, help="Output ACGT nucleotides instead of 0/1/2.")
    parser.add_argument("vcf", nargs="+", help="One or more VCFs (optionally gzipped)")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    genotypes = collections.defaultdict(list)
    all_samples = []

    def gt_map(s, alleles):
        # split on second char, e.g. "/" in "0/1/4"
        s_fields = s.split(s[1])
        if "." in s_fields:
            return "?"

        if args.nucleotide:
            # Use a random genotype.
            i = int(random.choice(s_fields))
            return alleles[i]

        alt_count = sum(1 for ss in s_fields if ss!="0")
        return str(alt_count)

    ref_list = []

    for filename in args.vcf:
        vparser = parsers.vcf_parser(filename)
        samples = next(vparser)
        all_samples.extend(samples)
        for row in vparser:
            (lineno, line, chrom, pos, id, ref, alt, qual, filter, info_str, fmt_fields, gts) = row

            alleles = [ref,]
            alleles.extend(alt.split(","))

            # Ignore indels, where REF or an ALT have length != 1.
            if len(alleles) != sum(len(a) for a in alleles):
                continue

            ref_list.append(ref)

            for sample in samples:
                gt = gt_map(gts[sample]["GT"], alleles)
                genotypes[sample].append(gt)

    nchars = len(ref_list)

    for sample in all_samples:
        if len(genotypes[sample]) != nchars:
            sys.stderr.write("Sample ``{}'' has {} characters, but reference has {}.\n".format(sample, len(genotypes[sample]), nchars))
            exit(1)

    if args.reference:
        ref_str = "hg19"
        all_samples.insert(0, ref_str)
        if args.nucleotide:
            genotypes[ref_str] = ref_list
        else:
            genotypes[ref_str] = "0"*len(ref_list)

    with sys.stdout as f:
        f.write("#NEXUS\n")
        f.write("BEGIN DATA;\n")
        f.write("DIMENSIONS NTAX={} NCHAR={};\n".format(len(all_samples), nchars))
        if args.nucleotide:
            f.write("FORMAT DATATYPE=DNA symbols=\"ACGT\" missing=? gap=- transpose=no;\n")
        else:
            f.write("FORMAT DATATYPE=integer symbols=\"012\" missing=? gap=- transpose=no;\n")
        f.write("matrix\n")

        for sample in all_samples:
            f.write("{} ".format(sample))
            f.write("".join(genotypes[sample]))
            f.write("\n")

        f.write(";\nend;\n")

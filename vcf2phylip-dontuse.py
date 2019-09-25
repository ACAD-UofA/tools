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
Convert one or more vcf to phylip format. SNPs are sorted by SNPid.

e.g.
$ ./vcf2phylip.py \
        --bed /localscratch/BisonProjects/SNPs/9908_Bovid_SNP_10k.bed \
        --min-depth 2 \
        --min-qual 25 \
        /localscratch/jsoubrier/10k/BisonX/run_raxml/*.vcf.gz \
        > Ancient_Bison.phy
"""

from __future__ import print_function
import sys
import collections
import parsers

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="Convert the SNPs from a vcf to phylip format (with 2=REF/REF, 1=REF/ALT, 0=ALT/ALT, ?=unknown).")
    parser.add_argument("--sortid", action="store_true", default=False, help="Sort bed file by SNP ID")
    parser.add_argument("--bed", required=True, help="Bed intervals to convert.")
    parser.add_argument("--min-depth", type=int, default=1, help="Filter calls based upon the DP field [%(default)s]")
    parser.add_argument("--min-qual", type=int, default=25, help="Filter calls based upon the QUAL field [%(default)s]")
    parser.add_argument("vcf", nargs="+", help="One or more VCFs (optionally gzipped)")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    # By default, we do *NOT* sort by SNP ID.  This is contrary to most of
    # the other scripts which deal in eigenstrat format.  The sort order
    # undertaken here is only to fit the sort order of existing data.
    loci_set, loci_list = parsers.parse_bed(args.bed, sort_by_id=args.sortid)

#    gt_map = {"0/0":"0", "0|0":"0",
#            "0/1":"1", "0|1":"1", "1/0":"1", "1|0":"1",
#            "1/1":"2", "1|1":"2",
#            "./.":"?", ".|.":"?"}

    gt_map = {"0/0":"2", "0|0":"2",
            "0/1":"1", "0|1":"1", "1/0":"1", "1|0":"1",
            "1/1":"0", "1|1":"0",
            "./.":"?", ".|.":"?"}
    genotypes = collections.defaultdict(dict)

    for filename in args.vcf:
        vparser = parsers.vcf_parser(filename)
        samples = next(vparser)
        for row in vparser:
            (lineno, line, chrom, pos, id, ref, alt, qual, filter, info_str, fmt_fields, gts) = row

            if (chrom, pos) not in loci_set:
                continue

            if float(qual) < args.min_qual:
                continue

            for sample in samples:
                dp = gts[sample].get("DP")
                if dp and int(dp) < args.min_depth:
                    continue
                gt = gt_map[gts[sample]["GT"]]
                if gt == "1":
                    continue
                genotypes[sample][(chrom, pos)] = gt

    sys.stdout.write(" {} {}\n".format(len(genotypes.keys()), len(loci_list)))

    for sample in genotypes.keys():
        sys.stdout.write("{} ".format(sample))
        sg_dict = genotypes[sample]
        gts = "".join(sg_dict.get((c, p), "?") for c, p, _ in loci_list)
        sys.stdout.write(gts)
        sys.stdout.write("\n")

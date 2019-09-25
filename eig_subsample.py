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
Subsample SNPs from eigensoft *.{snp,gno} files.

E.g.
$ eig_subsample.py \
        --snp 9908_749Bovid_5SteppeMN_1Caucasicus_5Steppe_9BisonX.snp \
        --geno 9908_749Bovid_5SteppeMN_1Caucasicus_5Steppe_9BisonX.geno \
        -f 0.5 \
        -o 50pc_subsample
"""

import random

def parse_lines(filename):
    with open(filename) as f:
        for line in f:
            yield line

def parse_args():
    import argparse

    parser = argparse.ArgumentParser(description="Subsample eigensoft *.{snp,geno} files.")
    parser.add_argument("--snp", required=True, help="eigensoft .snp file")
    parser.add_argument("--geno", required=True, help="eigensoft .geno file")
    parser.add_argument("-o", "--oprefix", required=True, help="Prefix for output *.{snp,geno} files.")
    parser.add_argument("-f", "--fraction", required=True, help="Fraction of the original snps to sample.")

    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    snps = list(parse_lines(args.snp))
    gts = list(parse_lines(args.geno))

    idx = random.sample(range(len(snps)), int(len(snps)*float(args.fraction)))

    with open("{}.snp".format(args.oprefix), "w") as f_snp, \
         open("{}.geno".format(args.oprefix), "w") as f_geno:
            for i in idx:
                f_snp.write(snps[i])
                f_geno.write(gts[i])

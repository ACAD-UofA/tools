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
Print stats on eigensoft .ind/.geno GT numbers.
"""

import sys
import re
import csv
import collections

def parse_ind(filename, accept_group=None):
    re_s = re.compile(r"[ \t]+")
    lineno = 0
    samples = []
    with open(filename) as f:
        for line in f:
            lineno += 1
            fields = re_s.split(line.strip())

            ind = fields[0]
            #sex = fields[1]
            group = fields[2]


            if accept_group and group != accept_group:
                continue

            samples.append((lineno-1, ind))

    return samples

def parse_geno(filename):
    with open(filename) as f:
        for line in f:
            yield line

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="Count genotypes for samples in eigensoft *.{ind,geno} files")
    parser.add_argument("--ind", required=True, help="eigensoft .ind file")
    parser.add_argument("--geno", required=True, help="eigensoft .geno file")
    parser.add_argument("--group", default=None, help="Limit output to individuals in the specified group")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    samples = parse_ind(args.ind, accept_group=args.group)
    count = collections.defaultdict(lambda:collections.Counter())

    for line in parse_geno(args.geno):
        for i, s in samples:
            gt = line[i]
            count[s][gt] += 1

    header = ["Sample", "Total SNPs", "REF/REF", "REF/ALT", "ALT/ALT"]
    writer = csv.DictWriter(sys.stdout, fieldnames=header, delimiter="\t")
    writer.writeheader()

    for _, s in samples:
        row = {}
        row["Sample"] = s
        row["Total SNPs"] = sum(count[s][k] for k in "012")
        row["REF/REF"] = count[s]["2"]
        row["REF/ALT"] = count[s]["1"]
        row["ALT/ALT"] = count[s]["0"]
        writer.writerow(row)

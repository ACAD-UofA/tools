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
Convert eigensoft *.{snp,ind,geno} files to Genepop format.
"""

import collections
import re
import os.path

def parse_lines(filename):
    with open(filename) as f:
        for line in f:
            yield line

def parse_ind(filename):
    re_s = re.compile(r"[ \t]+")
    groups = collections.defaultdict(list)
    for lineno, line in enumerate(parse_lines(filename)):
        fields = re_s.split(line.strip())

        ind = fields[0]
        #sex = fields[1]
        group = fields[2]

        groups[group].append((ind, lineno))

    return groups

def parse_snp(filename):
    re_s = re.compile(r"[ \t]+")
    snpids = []
    for lineno, line in enumerate(parse_lines(filename)):
        fields = re_s.split(line.strip())

        snpid = fields[0]
        #chrom = chrmap[fields[1]]
        #gpos = fields[2]
        #pos = int(fields[3])

        snpids.append(snpid)

    return snpids

def parse_args():
    import argparse

    parser = argparse.ArgumentParser(description="Convert eigensoft *.{snp,ind,geno} files to Genepop format.")
    parser.add_argument("--snp", required=True, help="eigensoft .snp file")
    parser.add_argument("--ind", required=True, help="eigensoft .ind file")
    parser.add_argument("--geno", required=True, help="eigensoft .geno file")
    parser.add_argument("--groups", help="Comma separated list of groups to include")

    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    # indicate the provenance of the data
    print("# eig2genepop.py - {}".format(os.path.basename(args.snp)[:-4]))

    snpids = parse_snp(args.snp)
    print(",".join(snpids))

    groups = parse_ind(args.ind)
    gtlines = list(parse_lines(args.geno))

    if args.groups:
        groupset = {g.strip() for g in args.groups.split(",")}

    gtmap = {"0": "0202", "1": "0102", "2": "0101", "9": "0000"}

    for group, indlist in groups.iteritems():
        if args.groups and group not in groupset:
            continue
        print("Pop")
        for ind, indno in indlist:
            gts = [gtmap[gtline[indno]] for gtline in gtlines]
            print("{}, {}".format(group, " ".join(gts)))
